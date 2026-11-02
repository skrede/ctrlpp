"""What each brace of a blanked C++ text opens, and the access sections and data members of every record."""

import bisect
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import NamedTuple

from style_source import spaced

DIRECTIVE = re.compile(r"^[ \t]*#(?:[^\n]*\\\n)*[^\n]*", re.M)
OPERATOR_SYMBOL = re.compile(r"\boperator\s*(\(\s*\)|\[\s*\]|<=>|<<=|>>=|->\*|<<|>>|->|&&|\|\||\+\+|--|[-+*/%^&|~!=<>,]=?)")
ATTRIBUTE = re.compile(r"\[\[[^\[\]]*\]\]")
LABELS = re.compile(r"(?:\s*(?:public|private|protected)\s*:(?!:))+")
LABEL = re.compile(r"\s*(public|private|protected)\s*:(?!:)")
TEMPLATE_CLAUSES = r"(?:template\s*<[^<>]*>\s*)*"
RECORD_HEAD = re.compile(r"\s*" + TEMPLATE_CLAUSES + r"(?:requires\b[^;{}]*?\s)?(class|struct|union)\b\s*"
                         r"(?:alignas\s*\([^()]*\)\s*)?(\w+(?:\s*::\s*\w+)*)?\s*(?:<[^<>]*>)?\s*(?:final\b)?\s*(?::(?!:).*)?",
                         re.S)
ENUM_HEAD = re.compile(r"\s*" + TEMPLATE_CLAUSES + r"(?:typedef\s+)?enum\b")
NAMESPACE_HEAD = re.compile(r'\s*(?:(?:inline\s+)?namespace\b|extern\s*"[^"]*"\s*$)')
GROUPED_SPECIFIER = re.compile(r"\b(?:decltype|alignas)\s*\(")
LAMBDA_SPECIFIER = re.compile(r"\s*(mutable|constexpr|consteval|static|noexcept)\b")
TOP_EQUAL = re.compile(r"(?<![=!<>])=(?!=)")
INIT_START = re.compile(r"(?<![=!<>])=(?!=)|\{")
SINGLE_COLON = re.compile(r"(?<!:):(?!:)")
COMPOUND_WORD = re.compile(r"\b(?:else|do|try)$")
REQUIRES_PREDECESSOR = re.compile(r"\b(?:requires|return)$")
POINTER_DECLARATOR = re.compile(r"\s*\*\s*(?:const\s+)?([A-Za-z_]\w*)\s*")
ARRAY_DIMENSIONS = re.compile(r"(?:\s*\[\s*\])+\s*$")
LAST_IDENTIFIER = re.compile(r"([A-Za-z_]\w*)\s*$")
IDENTIFIER_TAIL = re.compile(r"\w+$")
NOT_MEMBER_WORDS = frozenset({"using", "typedef", "friend", "static_assert", "template", "enum", "concept"})
RECORD_KEYWORDS = frozenset({"class", "struct", "union"})
LAMBDA_PREDECESSORS = frozenset({"return", "co_return", "co_yield", "throw"})
FUNCTION_CONTEXT = frozenset({"function", "block", "lambda"})
CONTINUING = frozenset({"initializer", "lambda"})
PAIRS = {"(": ")", "[": "]", "{": "}"}


@dataclass(eq=False)
class Block:
    kind: str
    open: int
    head: int
    parent: "Block | None"
    close: int = -1
    children: list["Block"] = field(default_factory=list)


class Scope(NamedTuple):
    blanked: str
    code: str
    root: Block
    newlines: list[int]

    def line(self, offset: int) -> int:
        return bisect.bisect_left(self.newlines, offset) + 1

    def innermost(self, offset: int) -> Block:
        block = self.root
        while True:
            index = bisect.bisect_left([child.open for child in block.children], offset) - 1
            if index < 0 or block.children[index].close <= offset:
                return block
            block = block.children[index]


class Statement(NamedTuple):
    start: int
    text: str
    access: str
    closer: Block | None


class Member(NamedTuple):
    name: str
    type: str
    access: str
    static: bool
    init: tuple[int, int] | None
    offset: int
    line: int


class Record(NamedTuple):
    kind: str
    name: str | None
    qualified: str
    nested: bool
    has_constructor: bool
    block: Block
    scope: Scope
    statements: tuple[Statement, ...]


def ident_before(text: str, index: int) -> bool:
    index -= 1
    while index >= 0 and text[index] in " \t":
        index -= 1
    return index >= 0 and (text[index].isalnum() or text[index] == "_")


def angle_close(text: str, start: int) -> int:
    depth, index = 0, start + 1
    while index < len(text):
        char = text[index]
        if char in "([":
            depth += 1
        elif char in ")]":
            if depth == 0:
                return -1
            depth -= 1
        elif char in ";{}":
            return -1
        elif char == "<" and ident_before(text, index):
            index = max(index, angle_close(text, index))
        elif char == ">" and depth == 0 and text[index - 1] != "-":
            return index
        index += 1
    return -1


def collapse_templates(text: str) -> str:
    """Blank every template argument list, keeping its angle brackets, so its commas and parentheses vanish."""
    pieces, copied, index = [], 0, 0
    while index < len(text):
        if text[index] == "<" and ident_before(text, index):
            close = angle_close(text, index)
            if close > 0:
                pieces += [text[copied:index + 1], spaced(text[index + 1:close])]
                copied, index = close, close
        index += 1
    return "".join(pieces) + text[copied:]


def masked(text: str, pattern: re.Pattern, group: int = 0, fill: str = " ") -> str:
    def replace(match: re.Match) -> str:
        whole, start = match.group(), match.start(group) - match.start()
        part = "".join(char if char == "\n" else fill for char in match.group(group))
        return whole[:start] + part + whole[start + len(part):]
    return pattern.sub(replace, text)


def prepare(blanked: str) -> str:
    return collapse_templates(masked(masked(blanked, DIRECTIVE), OPERATOR_SYMBOL, 1, "@"))


def matching(text: str, start: int) -> int:
    opening, closing, depth = text[start], PAIRS[text[start]], 0
    for index in range(start, len(text)):
        depth += (text[index] == opening) - (text[index] == closing)
        if depth == 0:
            return index
    return -1


def opening_of(text: str, close: int) -> int:
    depth = 0
    for index in range(close, -1, -1):
        depth += (text[index] == ")") - (text[index] == "(")
        if depth == 0:
            return index
    return -1


def shallow(text: str) -> str:
    """Blank the inside of every bracket group, keeping the outermost delimiters and every newline."""
    out, depth = [], 0
    for char in text:
        if char in "([{":
            out.append(char if depth == 0 else " ")
            depth += 1
        elif char in ")]}" and depth > 0:
            depth -= 1
            out.append(char if depth == 0 else " ")
        else:
            out.append(char if depth == 0 or char == "\n" else " ")
    return "".join(out)


def lambda_tail(text: str, index: int) -> bool:
    index = len(text) - len(text[index:].lstrip())
    if text.startswith("<", index):
        index = angle_close(text, index)
        if index < 0:
            return False
        index = len(text) - len(text[index + 1:].lstrip())
    if text.startswith("(", index):
        index = matching(text, index)
        if index < 0:
            return False
        index += 1
    while (specifier := LAMBDA_SPECIFIER.match(text, index)) is not None:
        index = len(text) - len(text[specifier.end():].lstrip())
        if specifier.group(1) == "noexcept" and text.startswith("(", index):
            index = matching(text, index) + 1
    rest = text[index:].strip()
    return not rest or rest.startswith("->") or rest.startswith("requires")


def is_lambda(head: str) -> bool:
    for start in (index for index, char in enumerate(head) if char == "[" and head[index + 1:index + 2] != "["):
        before = head[:start].rstrip()
        word = IDENTIFIER_TAIL.search(before)
        if (before and before[-1] in "])[>") or (word and word.group() not in LAMBDA_PREDECESSORS):
            continue
        close = matching(head, start)
        if close > 0 and lambda_tail(head, close + 1):
            return True
    return False


def is_requires_expression(head: str) -> bool:
    text = head.rstrip()
    if text.endswith(")"):
        opening = opening_of(text, len(text) - 1)
        if opening < 0:
            return False
        text = text[:opening].rstrip()
    if IDENTIFIER_TAIL.search(text) is None or IDENTIFIER_TAIL.search(text).group() != "requires":
        return False
    before = text[:-len("requires")].rstrip()
    return not before or before[-1] in "=&|(,!" or REQUIRES_PREDECESSOR.search(before) is not None


def is_compound(head: str) -> bool:
    text = head.strip()
    return not text or text.endswith(")") or (text.endswith(":") and not text.endswith("::")) or \
        COMPOUND_WORD.search(text) is not None


def declaration_kind(head: str) -> str:
    flat = shallow(head)
    paren, equal = flat.find("("), TOP_EQUAL.search(flat)
    if paren < 0 or (equal is not None and equal.start() < paren):
        return "initializer"
    if SINGLE_COLON.search(flat, flat.find(")", paren)) is None:
        return "function"
    tail = flat.rstrip()[-1:]
    return "initializer" if tail.isalnum() or tail in "_>" else "function"


def classify(head: str, parent: str) -> str:
    if parent == "requires":
        return "requires"
    if is_lambda(head):
        return "lambda"
    if parent in ("initializer", "enum"):
        return "initializer"
    if is_requires_expression(head):
        return "requires"
    if head.count("(") > head.count(")") or head.count("[") > head.count("]"):
        return "initializer"
    head = masked(head, ATTRIBUTE)
    head = head[LABELS.match(head).end():] if LABELS.match(head) else head
    if ENUM_HEAD.match(head):
        return "enum"
    if RECORD_HEAD.fullmatch(head):
        return "record"
    if parent in FUNCTION_CONTEXT:
        return "block" if is_compound(head) else "initializer"
    if NAMESPACE_HEAD.match(head):
        return "namespace"
    return declaration_kind(head)


@lru_cache(maxsize=4)
def parse(blanked: str) -> Scope:
    code = prepare(blanked)
    root = Block("file", -1, 0, None, len(code))
    stack, head = [root], 0
    for match in re.finditer(r"[;{}]", code):
        index, char = match.start(), match.group()
        if char == ";":
            head = index + 1
        elif char == "{":
            block = Block(classify(code[head:index], stack[-1].kind), index, head, stack[-1])
            stack[-1].children.append(block)
            stack.append(block)
            head = index + 1
        elif len(stack) > 1:
            block = stack.pop()
            block.close = index
            head = block.head if block.kind in CONTINUING else index + 1
    for block in stack[1:]:
        block.close = len(code)
    return Scope(blanked, code, root, [match.start() for match in re.finditer("\n", blanked)])


def walk(block: Block):
    for child in block.children:
        yield child
        yield from walk(child)


def flattened(code: str, block: Block) -> str:
    pieces, copied = [], block.open + 1
    for child in block.children:
        pieces += [code[copied:child.open + 1], spaced(code[child.open + 1:child.close])]
        copied = child.close
    return "".join(pieces) + code[copied:block.close]


def statements(scope: Scope, block: Block, access: str) -> tuple[Statement, ...]:
    """Split a record body at ';' and at each closing brace that ends a definition, tracking the access in force."""
    base, flat = block.open + 1, flattened(scope.code, block)
    closers = {child.close - base: child for child in block.children if child.kind not in CONTINUING}
    found, start = [], 0
    for index, char in enumerate(flat):
        if char == ";" or (char == "}" and index in closers):
            while (label := LABEL.match(flat, start)) is not None:
                access, start = label.group(1), label.end()
            stop = index if char == ";" else index + 1
            found.append(Statement(base + start, flat[start:stop], access, closers.get(index) if char == "}" else None))
            start = index + 1
    return tuple(found)


def cleaned(text: str) -> str:
    text = masked(text, ATTRIBUTE)
    for match in reversed(list(GROUPED_SPECIFIER.finditer(text))):
        close = matching(text, match.end() - 1)
        if close > 0:
            text = text[:match.end() - 1] + spaced(text[match.end() - 1:close + 1]) + text[close + 1:]
    return text


def pointer_declarator(text: str, flat: str) -> re.Match | None:
    paren = flat.find("(")
    if paren < 0:
        return None
    close = flat.find(")", paren)
    inner = POINTER_DECLARATOR.fullmatch(text, paren + 1, close)
    return inner if inner and flat[close + 1:].lstrip().startswith("(") else None


def is_function(text: str, closer: Block | None) -> bool:
    if closer is not None:
        return closer.kind == "function"
    text = cleaned(text)
    flat = shallow(text)
    paren, cut = flat.find("("), INIT_START.search(flat)
    return paren >= 0 and (cut is None or paren < cut.start()) and pointer_declarator(text, flat) is None


def declares_constructor(statement: Statement, name: str | None) -> bool:
    pattern = re.compile(r"(?<![\w~])(?<!operator )" + re.escape(name or "") + r"\s*\(")
    return name is not None and is_function(statement.text, statement.closer) and \
        pattern.search(shallow(cleaned(statement.text))) is not None


def record_of(scope: Scope, block: Block, qualifiers: dict, inherited: dict) -> Record:
    head = masked(scope.code[block.head:block.open], ATTRIBUTE)
    head = head[LABELS.match(head).end():] if LABELS.match(head) else head
    match = RECORD_HEAD.fullmatch(head)
    kind, name = match.group(1), match.group(2).split("::")[-1].strip() if match.group(2) else None
    nested = block.parent is not None and block.parent.kind == "record"
    access = "private" if kind == "class" else "public"
    access = inherited.get(block.open, access) if name is None and nested else access
    qualified = (qualifiers[block.parent] + "::" if nested else "") + (name or "<anonymous>")
    qualifiers[block] = qualified
    found = statements(scope, block, access)
    for statement in found:
        if statement.closer is not None:
            inherited[statement.closer.open] = statement.access
    return Record(kind, name, qualified, nested,
                  any(declares_constructor(statement, name) for statement in found), block, scope, found)


@lru_cache(maxsize=4)
def records(blanked: str) -> tuple[Record, ...]:
    scope, qualifiers, inherited = parse(blanked), {}, {}
    return tuple(record_of(scope, block, qualifiers, inherited) for block in walk(scope.root) if block.kind == "record")


def split_top_level(flat: str) -> list[tuple[int, int]]:
    cuts = [index for index, char in enumerate(flat) if char == ","]
    return list(zip([0] + [cut + 1 for cut in cuts], cuts + [len(flat)]))


def declarator(piece: tuple[int, int], text: str, flat: str):
    """Return (name offset, name, initializer span) for one declarator of a member declaration, or None."""
    begin, end = piece
    init = INIT_START.search(flat, begin, end)
    stop = init.start() if init else end
    pointer = pointer_declarator(text[:stop], flat[:stop]) if begin == 0 else None
    if pointer is None and SINGLE_COLON.search(flat, begin, stop):
        return None
    if pointer is not None:
        return pointer.start(1), pointer.group(1), (init.start(), end) if init else None
    declared = ARRAY_DIMENSIONS.sub("", flat[begin:stop])
    name = LAST_IDENTIFIER.search(declared)
    if name is None:
        return None
    return begin + name.start(1), name.group(1), (init.start(), end) if init else None


def statement_members(scope: Scope, statement: Statement) -> list[Member]:
    text = cleaned(statement.text)
    words = text.split()
    if statement.closer is not None or not words or words[0] in NOT_MEMBER_WORDS or \
            (words[0] in RECORD_KEYWORDS and len(words) == 2) or is_function(statement.text, None):
        return []
    flat, found, kind = shallow(text), [], None
    for piece in split_top_level(flat):
        parsed = declarator(piece, text, flat)
        if parsed is None:
            continue
        offset, name, init = parsed
        if kind is None:
            pointer = "(" in flat[:offset]
            kind = scope.blanked[statement.start:statement.start + ((init[0] if init else piece[1]) if pointer else offset)]
        if re.search(r"[A-Za-z_]", kind) is None:
            return []
        absolute = statement.start + offset
        span = (statement.start + init[0], statement.start + init[1]) if init else None
        found.append(Member(name, kind.strip(), statement.access, "static" in text[:offset].split(), span, absolute,
                            scope.line(absolute)))
    return found


def data_members(record: Record) -> list[Member]:
    return [member for statement in record.statements for member in statement_members(record.scope, statement)]
