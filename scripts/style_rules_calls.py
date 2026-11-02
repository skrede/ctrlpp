"""The discarded-result rule: a call to a fallible function made as a bare statement drops its typed failure."""

import re

from style_source import Source
from style_config import EXCLUDED_CALL_NAMES, FALLIBLE_CALL_NAMES
from style_scopes import FUNCTION_CONTEXT, Block, Scope, flattened, matching, parse, walk

Finding = tuple[int, str]

CALL_NAMES = tuple(name for name in FALLIBLE_CALL_NAMES if name not in EXCLUDED_CALL_NAMES)
BARE_CALL = re.compile(r"\s*(?:::\s*)?(?:[A-Za-z_]\w*\s*(?:<[^<>;{}]*>)?\s*(?:\.|->|::)\s*(?:template\s+)?)*"
                       r"(?P<name>" + "|".join(CALL_NAMES) + r")\b\s*(?:<[^<>;{}]*>)?\s*\(")
CONTROL_PREFIX = re.compile(r"\s*(?:(?:else|do|try)\b|case\b(?:[^;:]|::)*:(?!:)|default\s*:|"
                            r"(?P<open>(?:if|for|while|switch)\b\s*(?:constexpr\b)?\s*\())")
STATEMENT_BLOCKS = frozenset({"block", "function", "lambda", "record", "enum"})


def statements(scope: Scope, block: Block):
    """Yield (offset, text) for each statement of a function-context block that a ';' ends."""
    base, flat = block.open + 1, flattened(scope.code, block)
    edges = {edge - base for child in block.children if child.kind in STATEMENT_BLOCKS
             for edge in (child.open, child.close)}
    depth, start = 0, 0
    for index, char in enumerate(flat):
        depth += (char == "(") - (char == ")")
        if depth == 0 and char == ";":
            yield base + start, flat[start:index]
        if depth == 0 and (char == ";" or index in edges):
            start = index + 1


def after_control(text: str) -> int:
    index = 0
    while (prefix := CONTROL_PREFIX.match(text, index)) is not None:
        index = matching(text, prefix.end() - 1) + 1 if prefix.group("open") else prefix.end()
        if index <= 0:
            return len(text)
    return index


def discarded(text: str) -> re.Match | None:
    call = BARE_CALL.match(text, after_control(text))
    if call is None:
        return None
    close = matching(text, call.end() - 1)
    return call if close > 0 and not text[close + 1:].strip() else None


def check_discarded_result(source: Source) -> list[Finding]:
    scope, findings = parse(source.blanked), []
    for block in walk(scope.root):
        if block.kind not in FUNCTION_CONTEXT:
            continue
        for offset, text in statements(scope, block):
            if (call := discarded(text)) is not None:
                findings.append((scope.line(offset + call.start("name")),
                                 f"the result of {call.group('name')}(...) is discarded; bind, return or test it"))
    return findings
