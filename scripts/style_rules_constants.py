"""The constant-naming rule: namespace-scope constants, static constant members and constexpr locals are
lowercase with underscores, so a Hungarian k prefix or a camel-case hump fails."""

import re

from style_source import Source
from style_config import CONSTANT_CARVE_OUTS
from style_scopes import FUNCTION_CONTEXT, Block, Member, Scope, collapse_templates, parse, shallow, statement_members, \
    statements, walk

Finding = tuple[int, str]

HUNGARIAN = re.compile(r"k[A-Z]")
HUMP = re.compile(r"[a-z][A-Z]")
TEMPLATE_PREFIX = re.compile(r"\s*(?:template\s*<[^<>]*>\s*)+")
NAMESPACE_NAME = re.compile(r"namespace\s+(\w+(?:\s*::\s*\w+)*)")
RECORD_NAME = re.compile(r"\b(?:class|struct|union)\b\s*(?:alignas\s*\([^()]*\)\s*)?(\w+)")
CALLEE_NAME = re.compile(r"(~?\w+)\s*(?:<[^<>]*>)?\s*$")
NAMESPACE_KINDS = frozenset({"file", "namespace"})
CARVE_OUTS = tuple((prefix, re.compile(pattern)) for prefix, pattern, _ in CONSTANT_CARVE_OUTS)


def const_object(declared: str) -> bool:
    words = re.findall(r"\w+|[*&]", collapse_templates(declared))
    if "const" not in words:
        return False
    last_indirection = max((index for index, word in enumerate(words) if word in "*&"), default=-1)
    return max(index for index, word in enumerate(words) if word == "const") > last_indirection


def is_constant(member: Member, kind: str) -> bool:
    specifiers = member.type.split()
    if kind in NAMESPACE_KINDS:
        return "constexpr" in specifiers or const_object(member.type)
    if kind == "record":
        return member.static and ("constexpr" in specifiers or const_object(member.type))
    return "constexpr" in specifiers


def problem(name: str) -> str | None:
    if HUNGARIAN.match(name):
        return "a Hungarian k prefix"
    if HUMP.search(name):
        return "a camel-case hump"
    return None


def scope_name(scope: Scope, block: Block) -> str:
    names = []
    while block.parent is not None:
        head = shallow(scope.code[block.head:block.open])
        pattern = {"namespace": NAMESPACE_NAME, "record": RECORD_NAME}.get(block.kind)
        if pattern is not None:
            found = pattern.search(head)
            names.append(re.sub(r"\s+", "", found.group(1)) if found else "<anonymous>")
        elif block.kind == "function" and "(" in head and (found := CALLEE_NAME.search(head[:head.find("(")])):
            names.append(found.group(1))
        block = block.parent
    return "::".join(reversed(names))


def unprefixed(statement):
    prefix = TEMPLATE_PREFIX.match(statement.text)
    return statement._replace(text=" " * prefix.end() + statement.text[prefix.end():]) if prefix else statement


def constants(scope: Scope):
    for block in [scope.root, *walk(scope.root)]:
        if block.kind in NAMESPACE_KINDS or block.kind == "record" or block.kind in FUNCTION_CONTEXT:
            for statement in statements(scope, block, "public"):
                for member in statement_members(scope, unprefixed(statement)):
                    if is_constant(member, block.kind):
                        yield block, member


def carved_out(path: str, name: str) -> bool:
    return any(path.startswith(prefix) and pattern.fullmatch(name) for prefix, pattern in CARVE_OUTS)


def check_constant_naming(source: Source) -> list[Finding]:
    scope, findings = parse(source.blanked), []
    for block, member in constants(scope):
        if (why := problem(member.name)) is not None and not carved_out(source.path, member.name):
            qualified = "::".join(part for part in (scope_name(scope, block), member.name) if part)
            findings.append((member.line, f"{qualified}: a constant is lowercase with underscores ({why})"))
    return findings
