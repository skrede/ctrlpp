#!/usr/bin/env python3
"""Check every function in the tracked lib/ headers against the function-size ceiling and its register, as
clang-tidy's readability-function-size measures them."""

import argparse
import os
import shlex
import sys
from pathlib import Path

from check_file_size import PRECEDENCE, load_registry, report
from size_rules import REGISTRY_NAME, Key, Register, function_verdicts, is_cmake
from size_tidy import THRESHOLD, Finding, measure

EXTRA_ARG_VARIABLE = "CTRLPP_TIDY_EXTRA_ARG"

EXIT_HELP = f"""exit codes:
  0  every function over {THRESHOLD} body lines is registered, and every row still describes it
  1  a function is over the ceiling and not registered -- decompose it, or register it with a reason
  2  the arguments are wrong (argparse)
  3  a row is stale -- its function is at or under the ceiling or no longer exists; delete the row
  4  a registered function has grown past its {REGISTRY_NAME} allowance
  5  the register is missing or malformed
  6  clang-tidy failed, could not parse the headers, or did not report the generated canary function
  7  a row records more lines than its function has -- re-record it down

{EXTRA_ARG_VARIABLE}, when set, is split shell-style into compiler arguments passed before --extra-arg.
"""


def keyed(findings: list[Finding]) -> dict[Key, list[int]]:
    grouped = {}
    for path, name, lines in findings:
        grouped.setdefault((path, name), []).append(lines)
    return grouped


def cpp_rows(register: Register) -> dict[Key, list[int]]:
    return {key: rows for key, rows in register.functions.items() if not is_cmake(key[0])}


def run(args: argparse.Namespace) -> int:
    findings, source, failure = measure(args)
    if failure:
        return failure
    register, failure = load_registry(args.registry or args.source_root / REGISTRY_NAME, {})
    if failure:
        return failure
    rows = cpp_rows(register)
    result = function_verdicts(keyed(findings), rows)
    for code in PRECEDENCE:
        report(result[code])
    status = next((code for code in PRECEDENCE if result[code]), 0)
    if status == 0:
        print(f"checked {source}: {len(findings)} functions over {THRESHOLD} lines; "
              f"{sum(len(counts) for counts in rows.values())} registered functions within allowance")
    return status


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, epilog=EXIT_HELP,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add = parser.add_argument
    add("--source-root", type=Path, default=Path("."), help="repository root (default: the current directory)")
    add("--clang-tidy", default=os.environ.get("CLANG_TIDY") or "clang-tidy-18",
        help="clang-tidy binary (default: $CLANG_TIDY, else clang-tidy-18)")
    add("--eigen-include", type=Path, default=Path("/usr/include/eigen3"), help="Eigen include directory")
    add("--extra-arg", action="append", default=[], help="extra compiler argument (repeatable)")
    add("--registry", type=Path, default=None, help=f"register (default: {REGISTRY_NAME} under the root)")
    add("--tidy-output", type=Path, default=None, help="read a saved clang-tidy output instead of running it")
    args = parser.parse_args(argv)
    args.extra_arg = shlex.split(os.environ.get(EXTRA_ARG_VARIABLE, "")) + args.extra_arg
    return args


if __name__ == "__main__":
    raise SystemExit(run(parse_args(sys.argv[1:])))
