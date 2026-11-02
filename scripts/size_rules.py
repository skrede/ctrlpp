#!/usr/bin/env python3
"""The rules EXCEPTIONS.md encodes: which files it governs, what makes a row valid, and the verdicts a
measured unit earns against it."""

import re
from dataclasses import dataclass, field
from pathlib import Path

CEILING = 200
FUNCTION_CEILING = 25
GROWTH_NUMERATOR = 11
GROWTH_DENOMINATOR = 10
REGISTRY_NAME = "EXCEPTIONS.md"
CPP_PREFIX = "lib/"
GUIDELINE_PREFIXES = ("tests/", "examples/", "benchmarks/")
EXTENSIONS = frozenset({".h", ".hh", ".hpp", ".hxx", ".inl", ".ipp", ".c", ".cc", ".cpp", ".cxx"})
FILES_SECTION = "Files"
FUNCTIONS_SECTION = "Functions"

EXIT_UNLISTED = 1
EXIT_STALE = 3
EXIT_GROWTH = 4
EXIT_REGISTRY = 5
EXIT_NO_INPUT = 6
EXIT_INFLATED = 7

HEADING = re.compile(r"^#{1,6}\s+(.*?)\s*$")
SEPARATOR = re.compile(r"^\|(?:\s*:?-+:?\s*\|)+$")
COUNT_IN_REASON = re.compile(r"\b\d+\s*(?:lines?|-line)\b")
SPLIT_HARM = re.compile(r"split|decompos|separat|scatter", re.IGNORECASE)
FUTURE_WORK = re.compile(r"\b(?:planned|will be|to be split|todo|for now|eventually|later)\b",
                         re.IGNORECASE)
CMAKE_OPEN = re.compile(r"^\s*(function|macro)\s*\(\s*([^\s)]+)", re.IGNORECASE)
CMAKE_CLOSE = re.compile(r"^\s*end(function|macro)\s*\(", re.IGNORECASE)

Key = tuple[str, str]


@dataclass
class Register:
    files: dict[str, int] = field(default_factory=dict)
    functions: dict[Key, list[int]] = field(default_factory=dict)
    faults: list[str] = field(default_factory=list)


def is_cmake(path: str) -> bool:
    name = Path(path).name
    return name == "CMakeLists.txt" or name.endswith(".cmake") or name.endswith(".cmake.in")


def in_scope(path: str) -> bool:
    if is_cmake(path):
        return not path.startswith(GUIDELINE_PREFIXES)
    return path.startswith(CPP_PREFIX) and Path(path).suffix in EXTENSIONS


def allowance(recorded: int) -> int:
    return (recorded * GROWTH_NUMERATOR + GROWTH_DENOMINATOR - 1) // GROWTH_DENOMINATOR


def pair(measured: list[int], recorded: list[int]) -> tuple[list[tuple[int, int]], list[int], list[int]]:
    """Pair two multisets largest to largest; overloads share a key, so no row carries a line number."""
    ours, theirs = sorted(measured, reverse=True), sorted(recorded, reverse=True)
    n = min(len(ours), len(theirs))
    return list(zip(ours[:n], theirs[:n])), ours[n:], theirs[n:]


def cmake_function_lengths(text: str) -> list[tuple[str, int, int]]:
    """Each block spans its opening line to its closing line, the convention clang-tidy uses for braces."""
    found, open_blocks = [], []
    for number, line in enumerate(text.splitlines(), start=1):
        opening = CMAKE_OPEN.match(line)
        if opening:
            open_blocks.append((opening.group(2), number))
        elif CMAKE_CLOSE.match(line) and open_blocks:
            name, start = open_blocks.pop()
            found.append((name, start, number - start))
    return sorted(found, key=lambda block: block[1])


def function_verdicts(measured: dict[Key, list[int]], recorded: dict[Key, list[int]]) -> dict[int, list[str]]:
    """`measured` holds only the lengths past the function ceiling."""
    found = {EXIT_UNLISTED: [], EXIT_STALE: [], EXIT_GROWTH: [], EXIT_INFLATED: []}
    over = f"over the {FUNCTION_CEILING}-line ceiling"
    for key in sorted(set(measured) | set(recorded)):
        pairs, unlisted, stale = pair(measured.get(key, []), recorded.get(key, []))
        unit = f"{key[0]}:{key[1]}"
        beyond = "an overload beyond its rows in" if key in recorded else "and is not registered in"
        found[EXIT_UNLISTED] += [f"{unit}: {n} lines, {over}, {beyond} {REGISTRY_NAME}; decompose it"
                                 for n in unlisted]
        found[EXIT_STALE] += [f"{unit}: the {REGISTRY_NAME} row recording {n} lines is stale: the function "
                              f"is at or under the ceiling or gone; delete the row" for n in stale]
        found[EXIT_GROWTH] += [f"{unit}: {m} lines, {over}, grown from the {r} recorded in {REGISTRY_NAME} "
                               f"past its allowance of {allowance(r)}" for m, r in pairs if m > allowance(r)]
        found[EXIT_INFLATED] += [f"{unit}: the {REGISTRY_NAME} row records {r} lines but the function is "
                                 f"{m}; re-record it down" for m, r in pairs if r > m]
    return found


def registry_rows(lines: list[str]):
    """Every table line that is neither a separator nor a header is a row and must parse."""
    section = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        heading = HEADING.match(stripped)
        if heading:
            section = heading.group(1)
            continue
        if not stripped.startswith("|") or SEPARATOR.match(stripped):
            continue
        following = lines[index + 1].strip() if index + 1 < len(lines) else ""
        if SEPARATOR.match(following):
            continue
        cells = [cell.strip().strip("`").strip() for cell in stripped.strip("|").split("|")]
        yield index + 1, section, cells


def reason_fault(reason: str, counts: set[str], where: str) -> str | None:
    if not reason:
        return f"{where}: the reason is empty, which is a silent exception with extra steps"
    if COUNT_IN_REASON.search(reason) or any(re.search(rf"\b{n}\b", reason) for n in counts):
        return f"{where}: the reason restates the unit's line count instead of justifying it"
    if not SPLIT_HARM.search(reason):
        return f"{where}: the reason must name what splitting this unit would harm"
    if FUTURE_WORK.search(reason):
        return f"{where}: the reason promises future work; this registers decisions already taken"
    return None


def file_row_fault(cells: list[str], where: str, counts: dict[str, int]) -> str | None:
    if len(cells) < 3 or not cells[0]:
        return f"{where}: a file row needs a path, a recorded line count and a reason"
    if not in_scope(cells[0]):
        return f"{where}: the path is outside the checked scope"
    if not cells[1].isdigit():
        return f"{where}: the recorded count column is not a number"
    measured = counts.get(cells[0])
    known = {cells[1]} | ({str(measured)} if measured is not None else set())
    return reason_fault(cells[2], known, where)


def function_row_fault(cells: list[str], where: str) -> str | None:
    if len(cells) < 4 or not cells[0] or not cells[1]:
        return f"{where}: a function row needs a path, a function, a recorded line count and a reason"
    if not in_scope(cells[0]):
        return f"{where}: the path is outside the checked scope"
    if not cells[2].isdigit():
        return f"{where}: the recorded count column is not a number"
    return reason_fault(cells[3], {cells[2]}, where)


def add_file_row(register: Register, cells: list[str], where: str, counts: dict[str, int]) -> None:
    fault = file_row_fault(cells, where, counts)
    if fault is not None:
        register.faults.append(fault)
    elif cells[0] in register.files:
        register.faults.append(f"{where}: the path is registered twice")
    else:
        register.files[cells[0]] = int(cells[1])


def add_function_row(register: Register, cells: list[str], where: str) -> None:
    fault = function_row_fault(cells, where)
    if fault is not None:
        register.faults.append(fault)
    else:
        register.functions.setdefault((cells[0], cells[1]), []).append(int(cells[2]))


def parse_registry(text: str, counts: dict[str, int]) -> Register:
    register = Register()
    for number, section, cells in registry_rows(text.splitlines()):
        where = f"{REGISTRY_NAME} line {number}, {cells[0] or '<no path>'}"
        if section == FILES_SECTION:
            add_file_row(register, cells, where, counts)
        elif section == FUNCTIONS_SECTION:
            add_function_row(register, cells, where)
        else:
            register.faults.append(f"{where}: a table row under the heading {section!r}, which is "
                                   f"neither {FILES_SECTION!r} nor {FUNCTIONS_SECTION!r}")
    return register
