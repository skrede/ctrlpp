"""Tracked-file enumeration and offset-preserving blanking of C++ comments and literals for the style lint."""

import re
import subprocess
from pathlib import Path
from typing import NamedTuple

from style_config import CPP_SCOPE_PREFIXES, CPP_SUFFIXES, SPELLING_SUFFIXES, SPELLING_NAMES

RAW_PREFIXES = ("R", "u8R", "uR", "UR", "LR")
RAW_OPENING = re.compile(r'"([^()\\\s"]{0,16})\(')


class Source(NamedTuple):
    path: str
    raw: str
    blanked: str
    tracked: frozenset[str]


def tracked_files(root: Path) -> list[str]:
    listing = subprocess.run(["git", "-C", str(root), "ls-files", "-z"], check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    return sorted(name for name in listing.decode("utf-8", "surrogateescape").split("\0") if name)


def in_cpp_scope(path: str) -> bool:
    return path.startswith(CPP_SCOPE_PREFIXES) and path.endswith(CPP_SUFFIXES)


def in_spelling_scope(path: str) -> bool:
    return in_cpp_scope(path) or path.rsplit("/", 1)[-1] in SPELLING_NAMES or path.endswith(SPELLING_SUFFIXES)


def line_of(text: str, index: int) -> int:
    return text.count("\n", 0, index) + 1


def spaced(text: str) -> str:
    return "".join(char if char == "\n" else " " for char in text)


def end_of_line_comment(text: str, start: int) -> int:
    index = start
    while True:
        index = text.find("\n", index)
        if index < 0:
            return len(text)
        if text[index - 1] != "\\":
            return index
        index += 1


def end_of_quoted(text: str, start: int) -> int:
    quote, index = text[start], start + 1
    while index < len(text) and text[index] != quote and text[index] != "\n":
        index += 2 if text[index] == "\\" else 1
    return min(index, len(text))


def token_before(text: str, index: int, extra: str = "") -> str:
    start = index
    while start > 0 and (text[start - 1].isalnum() or text[start - 1] in "_" + extra):
        start -= 1
    return text[start:index]


def is_digit_separator(text: str, index: int) -> bool:
    token = token_before(text, index, ".'")
    return bool(token) and token[0].isdigit() and index + 1 < len(text) and text[index + 1].isalnum()


def raw_string_close(text: str, quote: int) -> int:
    if token_before(text, quote) not in RAW_PREFIXES:
        return -1
    opening = RAW_OPENING.match(text, quote)
    if not opening:
        return -1
    closing = text.find(")" + opening.group(1) + '"', opening.end())
    return len(text) if closing < 0 else closing + len(opening.group(1)) + 1


def closed(text: str, close: int) -> int:
    return close + 1 if close < len(text) and text[close] in "\"'" else close


def literal_span(text: str, index: int) -> tuple[int, int, int]:
    """Return (body start, body end, resume) for a comment or literal opening at index, else an empty body."""
    char, pair = text[index], text[index:index + 2]
    if pair == "//":
        end = end_of_line_comment(text, index)
        return index, end, end
    if pair == "/*":
        close = text.find("*/", index + 2)
        end = len(text) if close < 0 else close + 2
        return index, end, end
    if char == '"':
        close = raw_string_close(text, index)
        close = end_of_quoted(text, index) if close < 0 else close
        return index + 1, close, closed(text, close)
    if char == "'" and not is_digit_separator(text, index):
        close = end_of_quoted(text, index)
        return index + 1, close, closed(text, close)
    return index, index, index + 1


def blank(text: str) -> str:
    """Replace every comment and every string, raw-string and character literal body with spaces.

    The result has the length and newline positions of the input, so an index or a line number found in it
    holds in the original text. Literal delimiters stay, so an empty literal still reads as a token.
    """
    pieces, index, copied = [], 0, 0
    while index < len(text):
        start, end, resume = literal_span(text, index)
        if end > start:
            pieces += [text[copied:start], spaced(text[start:end])]
            copied = end
        index = resume
    pieces.append(text[copied:])
    return "".join(pieces)
