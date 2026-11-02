#!/usr/bin/env python3
"""Check every tracked lib/ C++ file and every tracked CMake file outside the test, example and benchmark
trees against the file-size ceiling, and every CMake function and macro against the function-size ceiling,
both against their register."""

import argparse
import subprocess
import sys
from pathlib import Path

from size_rules import (CEILING, EXIT_GROWTH, EXIT_INFLATED, EXIT_NO_INPUT, EXIT_REGISTRY, EXIT_STALE,
                        EXIT_UNLISTED, FUNCTION_CEILING, REGISTRY_NAME, Key, Register, allowance,
                        cmake_function_lengths, function_verdicts, in_scope, is_cmake, parse_registry)

PRECEDENCE = (EXIT_UNLISTED, EXIT_STALE, EXIT_GROWTH, EXIT_INFLATED)

EXIT_HELP = f"""exit codes:
  0  every file over {CEILING} lines and every CMake function over {FUNCTION_CEILING} is registered, and every
     row still describes its unit
  1  a unit is over its ceiling and not registered -- decompose it, or register it with a reason
  2  the arguments are wrong (argparse)
  3  a row is stale -- its unit is at or under the ceiling or no longer exists; delete the row
  4  a registered unit has grown past its {REGISTRY_NAME} allowance
  5  the register is missing or malformed
  6  the scan could not read what it was asked to read
  7  a row records more lines than its unit has -- re-record it down
"""


def tracked_sources(root: Path) -> list[str]:
    listing = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout
    names = listing.decode("utf-8", "surrogateescape").split("\0")
    return sorted(name for name in names if name and in_scope(name))


def read_text(path: Path) -> str:
    return path.read_bytes().decode("utf-8", "surrogateescape")


def line_count(text: str) -> int:
    return text.count("\n") + (1 if text and not text.endswith("\n") else 0)


def measure(root: Path, names: list[str]) -> tuple[dict[str, int], dict[Key, list[int]], list[str]]:
    counts, functions, absent = {}, {}, []
    for name in names:
        path = root / name
        if not path.is_file():
            absent.append(name)
            continue
        text = read_text(path)
        counts[name] = line_count(text)
        for function, _, length in cmake_function_lengths(text) if is_cmake(name) else []:
            if length > FUNCTION_CEILING:
                functions.setdefault((name, function), []).append(length)
    return counts, functions, absent


def report(messages) -> None:
    for message in messages:
        print(f"error: {message}", file=sys.stderr)


def scan(root: Path) -> tuple[dict[str, int], dict[Key, list[int]], int]:
    try:
        names = tracked_sources(root)
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"error: cannot list the tracked files under {root}: {exc}", file=sys.stderr)
        return {}, {}, EXIT_NO_INPUT
    counts, functions, absent = measure(root, names)
    report(f"{path} is tracked but absent from the working tree, so its length cannot be measured; "
           "restore it, or remove it from the index" for path in absent)
    if not counts:
        print(f"error: no tracked file in the checked scope was found under {root}", file=sys.stderr)
    return counts, functions, EXIT_NO_INPUT if absent or not counts else 0


def load_registry(path: Path, counts: dict[str, int]) -> tuple[Register, int]:
    if not path.is_file():
        print(f"error: the register {path} does not exist", file=sys.stderr)
        return Register(), EXIT_REGISTRY
    register = parse_registry(path.read_text(encoding="utf-8"), counts)
    report(register.faults)
    return register, EXIT_REGISTRY if register.faults else 0


def exception_findings(counts: dict[str, int], recorded: dict[str, int]) -> dict[int, list[str]]:
    live = {p: counts[p] for p in recorded if counts.get(p, 0) > CEILING}
    return {
        EXIT_STALE: [f"{p}: the {REGISTRY_NAME} row is stale: the file is at or under the {CEILING}-line "
                     "ceiling or gone; delete the row" for p in sorted(recorded) if p not in live],
        EXIT_GROWTH: [f"{p} has grown to {n} lines from the {recorded[p]} recorded in {REGISTRY_NAME}, "
                      f"past its allowance of {allowance(recorded[p])}"
                      for p, n in sorted(live.items()) if n > allowance(recorded[p])],
        EXIT_INFLATED: [f"{p}: the {REGISTRY_NAME} row records {recorded[p]} lines but the file has {n}; "
                        "re-record it down" for p, n in sorted(live.items()) if recorded[p] > n],
    }


def cmake_rows(register: Register) -> dict[Key, list[int]]:
    return {key: rows for key, rows in register.functions.items() if is_cmake(key[0])}


def verdicts(counts: dict[str, int], functions: dict[Key, list[int]], register: Register) -> dict[int, list[str]]:
    findings = {EXIT_UNLISTED: [f"{p} has {n} lines, over the {CEILING}-line ceiling, and is not registered "
                                f"in {REGISTRY_NAME}" for p, n in sorted(counts.items())
                                if n > CEILING and p not in register.files]}
    for code, messages in exception_findings(counts, register.files).items():
        findings[code] = messages
    for code, messages in function_verdicts(functions, cmake_rows(register)).items():
        findings[code] += messages
    return findings


def run(args: argparse.Namespace) -> int:
    counts, functions, failure = scan(args.source_root)
    if failure:
        return failure
    register, failure = load_registry(args.registry or args.source_root / REGISTRY_NAME, counts)
    if failure:
        return failure
    findings = verdicts(counts, functions, register)
    for code in PRECEDENCE:
        report(findings[code])
    status = next((code for code in PRECEDENCE if findings[code]), 0)
    if status == 0:
        cmake = sum(1 for name in counts if is_cmake(name))
        rows = sum(len(lengths) for lengths in cmake_rows(register).values())
        print(f"checked {len(counts)} tracked files ({cmake} CMake); {len(register.files)} registered files "
              f"and {rows} registered functions within allowance")
    return status


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, epilog=EXIT_HELP,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source-root", type=Path, default=Path("."),
                        help="repository root to check (default: the current directory)")
    parser.add_argument("--registry", type=Path, default=None,
                        help=f"register to read (default: {REGISTRY_NAME} under the source root)")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(run(parse_args(sys.argv[1:])))
