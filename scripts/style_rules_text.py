"""Style rules decided line by line: attributes, #pragma once, header guards, closing comments, spelling."""

import re

from style_source import Source, line_of
from style_config import (ALLOWED_ATTRIBUTES, ALLOWED_PHRASES, ATTRIBUTE_CARVE_OUTS, BRITISH_PATTERNS, GUARD_PREFIXES,
                          HEADER_SUFFIXES)

Finding = tuple[int, str]

# The closing ']]' must follow the opening '[[' with no bracket between them, which is what tells an attribute
# apart from an immediately invoked lambda used as a subscript, a[[&]{ return 0; }()], whose capture closes on a
# single ']'.
ATTRIBUTE_SYNTAX = re.compile(r"\[\[([^\[\]]*)\]\]")
ATTRIBUTE_USING = re.compile(r"^\s*using\s+(\w+)\s*:(.*)$", re.S)
ATTRIBUTE_NAME = re.compile(r"^\s*(\w+(?:::\w+)?)")
PRAGMA_ONCE = re.compile(r"^[ \t]*#[ \t]*pragma[ \t]+once\b", re.M)
DIRECTIVE = re.compile(r"^[ \t]*#[ \t]*(\w+)[ \t]*(\w*)", re.M)
NAMESPACE_CLOSER = re.compile(r"\}[ \t]*(?://|/\*)[ \t]*(?:end[ \t]+(?:of[ \t]+)?)?namespace\b", re.I)
URL = re.compile(r"\b[a-z][a-z0-9+.-]*://\S+|\bwww\.\S+", re.I)
BRITISH = tuple((re.compile(pattern, re.I), american) for pattern, american in BRITISH_PATTERNS)


def top_level_split(content: str) -> list[str]:
    parts, depth, start = [], 0, 0
    for index, char in enumerate(content):
        depth += (char == "(") - (char == ")")
        if char == "," and depth == 0:
            parts.append(content[start:index])
            start = index + 1
    return parts + [content[start:]]


def attribute_names(content: str) -> list[str]:
    using = ATTRIBUTE_USING.match(content)
    scope, body = (using.group(1) + "::", using.group(2)) if using else ("", content)
    names = (ATTRIBUTE_NAME.match(part) for part in top_level_split(body))
    return [scope + name.group(1) for name in names if name]


def attribute_allowed(path: str, name: str) -> bool:
    return name in ALLOWED_ATTRIBUTES or any(
        path.startswith(prefix) and name == carved for prefix, carved, _ in ATTRIBUTE_CARVE_OUTS)


def check_attributes(source: Source) -> list[Finding]:
    return [(line_of(source.blanked, match.start()),
             f"[[{name}]] is an attribute the compiler does not require; remove it")
            for match in ATTRIBUTE_SYNTAX.finditer(source.blanked)
            for name in attribute_names(match.group(1)) if not attribute_allowed(source.path, name)]


def check_pragma_once(source: Source) -> list[Finding]:
    return [(line_of(source.blanked, match.start()), "#pragma once is not used; guard the header instead")
            for match in PRAGMA_ONCE.finditer(source.blanked)]


def expected_guard(path: str) -> str | None:
    for prefix, guard in GUARD_PREFIXES:
        if path.startswith(prefix):
            return guard + re.sub(r"[^A-Za-z0-9]", "_", path[len(prefix):]).upper()
    return None


def guard_endif(blanked: str) -> re.Match | None:
    """The guard's #endif: the last directive, standing at the end of a header opened by #ifndef G / #define G."""
    directives = list(DIRECTIVE.finditer(blanked))
    if len(directives) < 3 or blanked[:directives[0].start()].strip():
        return None
    opening, definition, last = directives[0], directives[1], directives[-1]
    if (opening.group(1), definition.group(1), last.group(1)) != ("ifndef", "define", "endif"):
        return None
    if opening.group(2) != definition.group(2) or blanked[last.end():].strip():
        return None
    return last


def check_guard(source: Source) -> list[Finding]:
    guard = expected_guard(source.path)
    if guard is None or not source.path.endswith(HEADER_SUFFIXES):
        return []
    if guard_endif(source.blanked) is None:
        return [(1, f"the header is not guarded as a whole; open it with #ifndef {guard} and #define {guard} "
                    "and close it with #endif")]
    opening = DIRECTIVE.search(source.blanked)
    if opening.group(2) != guard:
        return [(line_of(source.blanked, opening.start()), f"the guard is {opening.group(2)}; it should be {guard}")]
    return []


def check_closing_comment(source: Source) -> list[Finding]:
    findings = [(line_of(source.raw, match.start()), "a comment after a closing namespace brace; remove it")
                for match in NAMESPACE_CLOSER.finditer(source.raw) if source.blanked[match.start()] == "}"]
    endif = guard_endif(source.blanked) if source.path.endswith(HEADER_SUFFIXES) else None
    if endif is not None:
        trailing = re.search(r"\S", source.raw[endif.end(1):])
        if trailing:
            findings.append((line_of(source.raw, endif.end(1) + trailing.start()),
                             "a comment after the guard's #endif; remove it"))
    return findings


def masked(text: str, span: tuple[int, int]) -> str:
    return text[:span[0]] + " " * (span[1] - span[0]) + text[span[1]:]


def spelling_text(raw: str) -> str:
    text = raw
    for match in URL.finditer(raw):
        text = masked(text, match.span())
    for phrase, _ in ALLOWED_PHRASES:
        for match in re.finditer(re.escape(phrase), raw):
            text = masked(text, match.span())
    return text


def check_spelling(source: Source) -> list[Finding]:
    text = spelling_text(source.raw)
    findings = [(line_of(text, match.start()), f"'{match.group()}' is British spelling; write '{match.expand(american)}'")
                for pattern, american in BRITISH for match in pattern.finditer(text)]
    return sorted(findings)
