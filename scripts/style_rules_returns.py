"""The return-type rule: a function declared at namespace or record scope spells its return type in front."""

import re

from style_source import Source
from style_scopes import matching, parse
from style_config import UNSPELLABLE_RETURNS

Finding = tuple[int, str]

AUTO_HEAD = re.compile(r"\b(?:decltype\s*\(\s*auto\s*\)|auto)(?:\s*(?:\bconst\b|\bvolatile\b|&&|&|\*))*\s*"
                       r"(?P<name>(?:[A-Za-z_]\w*\s*(?:<[^<>;{}]*>)?\s*::\s*)*(?:operator\s*@+|[A-Za-z_]\w*))\s*\(")
SPECIFIERS = frozenset({"static", "constexpr", "consteval", "constinit", "inline", "explicit", "friend", "virtual",
                        "extern", "const", "volatile"})
QUALIFIER = re.compile(r"\s*(const\b|volatile\b|&&|&|noexcept\b|throw\s*\(\s*\))")
TYPE_END = re.compile(r"[{;]|(?<![=!<>])=(?!=)|\b(?:requires|override|final)\b")
OPEN_PAREN = re.compile(r"\s*\(")
ARROW = re.compile(r"\s*->")
SCOPES = frozenset({"file", "namespace", "record"})


def starts_declaration(code: str, index: int) -> bool:
    """Whether only specifiers, attributes and whitespace separate a head from the end of what precedes it."""
    while True:
        end = index
        while end > 0 and code[end - 1].isspace():
            end -= 1
        word = end
        while word > 0 and (code[word - 1].isalnum() or code[word - 1] == "_"):
            word -= 1
        if word < end and code[word:end] in SPECIFIERS:
            index = word
        elif code[max(0, end - 2):end] == "]]" and code.rfind("[[", 0, end) >= 0:
            index = code.rfind("[[", 0, end)
        else:
            return end == 0 or (code[end - 1] in ";{}>):" and code[max(0, end - 2):end] != "->")


def trailing_type(code: str, blanked: str, close: int) -> str | None:
    index = close + 1
    while (qualifier := QUALIFIER.match(code, index)) is not None:
        index = qualifier.end()
        group = OPEN_PAREN.match(code, index)
        if qualifier.group(1) == "noexcept" and group is not None:
            index = matching(code, group.end() - 1) + 1
    arrow = ARROW.match(code, index)
    if arrow is None:
        return None
    end = TYPE_END.search(code, arrow.end())
    return " ".join(blanked[arrow.end():end.start() if end else len(code)].split())


def unspellable(path: str, name: str) -> bool:
    return any(path == listed and name == function for listed, function, _ in UNSPELLABLE_RETURNS)


def check_return_type(source: Source) -> list[Finding]:
    scope, findings = parse(source.blanked), []
    for head in AUTO_HEAD.finditer(scope.code):
        name = re.sub(r"\s+", "", source.blanked[head.start("name"):head.end("name")])
        if scope.innermost(head.start()).kind not in SCOPES or not starts_declaration(scope.code, head.start()) or \
                unspellable(source.path, name.split("::")[-1]):
            continue
        returned = trailing_type(scope.code, source.blanked, matching(scope.code, head.end() - 1))
        if returned is None:
            findings.append((scope.line(head.start()), f"{name}: deduced return type; spell it"))
        elif not re.match(r"decltype\s*\((?!\s*auto\s*\))", returned):
            findings.append((scope.line(head.start()), f"{name}: a trailing return type counts as auto; spell the "
                                                       f"return type {returned} in front"))
    return findings
