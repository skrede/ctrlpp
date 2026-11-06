#!/usr/bin/env python3
"""Lint the tracked C++ files under lib/ and examples/, and the tracked CMake files for spelling, against the
mechanical rules of CONVENTIONS.md. Comments and literals are blanked before every rule but spelling."""

import argparse
import posixpath
import subprocess
import sys
from pathlib import Path
from typing import Callable, NamedTuple

from style_rules_calls import check_discarded_result
from style_rules_returns import check_return_type
from style_rules_includes import check_include_order
from style_rules_constants import check_constant_naming
from style_source import Source, blank, in_cpp_scope, in_spelling_scope, tracked_files
from style_rules_members import check_member_init, check_member_naming, list_members
from style_rules_text import check_attributes, check_closing_comment, check_guard, check_pragma_once, check_spelling

EXIT_USAGE = 2
EXIT_NO_INPUT = 6


class Rule(NamedTuple):
    code: int
    check: Callable[[Source], list[tuple[int, str]]]
    applies: Callable[[str], bool]


RULES = {
    "include-order": Rule(10, check_include_order, in_cpp_scope),
    "attributes": Rule(11, check_attributes, in_cpp_scope),
    "return-type": Rule(12, check_return_type, in_cpp_scope),
    "member-init": Rule(13, check_member_init, in_cpp_scope),
    "member-naming": Rule(14, check_member_naming, in_cpp_scope),
    "pragma-once": Rule(15, check_pragma_once, in_cpp_scope),
    "guard": Rule(16, check_guard, in_cpp_scope),
    "closing-comment": Rule(17, check_closing_comment, in_cpp_scope),
    "spelling": Rule(18, check_spelling, in_spelling_scope),
    "discarded-result": Rule(19, check_discarded_result, in_cpp_scope),
    "constant-naming": Rule(20, check_constant_naming, in_cpp_scope),
}

EXIT_HELP = "exit codes:\n  0  no findings\n" + "".join(
    f"  {rule.code}  a {name} finding (the lowest code among the failing rules is returned)\n"
    for name, rule in RULES.items()) + (f"  {EXIT_USAGE}  the arguments are wrong, or a path is outside the lint's scope\n"
                                       f"  {EXIT_NO_INPUT}  a file could not be read, or there was nothing to check\n")


def read_text(path: Path) -> str:
    return path.read_bytes().decode("utf-8", "surrogateescape")


def requested(root: Path, paths: list[str], parser: argparse.ArgumentParser) -> list[str]:
    names = []
    for given in paths:
        path = Path(given)
        try:
            name = path.resolve().relative_to(root.resolve()).as_posix() if path.is_absolute() else \
                posixpath.normpath(path.as_posix())
        except ValueError:
            name = given
        if not in_spelling_scope(name):
            parser.error(f"{given} is outside the lint's scope (C++ files under lib/ and examples/, and CMake files)")
        names.append(name)
    return names


def lint(root: Path, names: list[str], tracked: frozenset[str], rules: dict[str, Rule]) -> tuple[list, int, int]:
    findings, checked = [], 0
    for name in names:
        selected = [(rule_id, rule) for rule_id, rule in rules.items() if rule.applies(name)]
        if not selected:
            continue
        try:
            raw = read_text(root / name)
        except OSError as exc:
            print(f"error: cannot read {name}: {exc.strerror}", file=sys.stderr)
            return findings, checked, EXIT_NO_INPUT
        source = Source(name, raw, blank(raw) if in_cpp_scope(name) else raw, tracked)
        checked += 1
        findings += [(name, line, rule.code, rule_id, message)
                     for rule_id, rule in selected for line, message in rule.check(source)]
    return findings, checked, 0


def summary(checked: int, findings: list) -> str:
    if not findings:
        return f"checked {checked} files, no findings"
    counts = {}
    for _, _, code, rule_id, _ in sorted(findings, key=lambda finding: finding[2]):
        counts[rule_id] = counts.get(rule_id, 0) + 1
    tally = ", ".join(f"{rule_id}: {count}" for rule_id, count in counts.items())
    return f"checked {checked} files; {len(findings)} findings ({tally})"


def print_members(root: Path, names: list[str]) -> int:
    for name in (name for name in names if in_cpp_scope(name)):
        try:
            raw = read_text(root / name)
        except OSError as exc:
            print(f"error: cannot read {name}: {exc.strerror}", file=sys.stderr)
            return EXIT_NO_INPUT
        print("".join(line + "\n" for line in list_members(Source(name, raw, blank(raw), frozenset()))), end="")
    return 0


def run(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.list_rules:
        print("".join(f"{name} {rule.code}\n" for name, rule in RULES.items()), end="")
        return 0
    if args.list_members and args.rule != ["member-init"]:
        parser.error("--list-members goes with --rule member-init and no other rule")
    rules = {name: RULES[name] for name in RULES if not args.rule or name in args.rule}
    names = requested(args.source_root, args.paths, parser)
    try:
        tracked = frozenset(tracked_files(args.source_root))
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"error: cannot list the tracked files under {args.source_root}: {exc}", file=sys.stderr)
        return EXIT_NO_INPUT
    if args.list_members:
        return print_members(args.source_root, names or sorted(tracked))
    findings, checked, failure = lint(args.source_root, names or sorted(tracked), tracked, rules)
    if failure:
        return failure
    if checked == 0:
        print(f"error: nothing to check under {args.source_root}", file=sys.stderr)
        return EXIT_NO_INPUT
    for path, line, _, rule_id, message in sorted(findings):
        print(f"error: {path}:{line}: {rule_id}: {message}", file=sys.stderr)
    print(summary(checked, findings))
    return min((finding[2] for finding in findings), default=0)


def parser_for() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, epilog=EXIT_HELP,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source-root", type=Path, default=Path("."),
                        help="repository root to lint (default: the current directory)")
    parser.add_argument("--rule", action="append", choices=list(RULES), metavar="ID",
                        help="run only this rule; repeat for more (default: every rule)")
    parser.add_argument("--list-rules", action="store_true", help="print each rule's id and exit code, then exit")
    parser.add_argument("--list-members", action="store_true",
                        help="with --rule member-init, print each member the rule flags with its initializer, "
                             "then exit 0")
    parser.add_argument("paths", nargs="*", metavar="PATH",
                        help="files to lint, relative to the source root (default: every tracked file in scope)")
    return parser


if __name__ == "__main__":
    PARSER = parser_for()
    raise SystemExit(run(PARSER.parse_args(), PARSER))
