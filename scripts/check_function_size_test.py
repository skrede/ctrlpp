#!/usr/bin/env python3
"""Fixture repositories and canned clang-tidy output exercising every verdict the function-size gate can reach."""

import subprocess
import sys
import tempfile
from pathlib import Path

GATE = Path(__file__).resolve().parent / "check_function_size.py"
WHY = "one cohesive routine whose shared locals splitting would scatter across helpers"
FILES = "# register\n\n## Files\n\n| Path | Lines at registration | Why |\n| --- | --- | --- |\n"
FUNCTIONS = "\n## Functions\n\n| Path | Function | Lines at registration | Why |\n| --- | --- | --- | --- |\n"
CHECK = "[readability-function-size]"
HEADER = "lib/ctrlpp/include/ctrlpp/control/dare.h"
OTHER = "lib/ctrlpp/include/ctrlpp/detail/kernel.h"
ERROR = "/x/tu.cpp:1:10: error: 'missing.h' file not found [clang-diagnostic-error]\n"
ORPHAN = f"/x/{HEADER}:9:5: warning: function 'orphan' exceeds recommended size/complexity thresholds {CHECK}\n"


def git(root: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(root), "-c", "user.name=fixture", "-c", "user.email=f@x",
                           "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, text=True).stdout


def diagnostic(path: str, line: int, name: str, lines: int) -> str:
    at = f"{path}:{line}:5"
    return (f"{at}: warning: function '{name}' exceeds recommended size/complexity thresholds {CHECK}\n"
            f"    void {name}()\n    ^\n"
            f"{at}: note: {lines} lines including whitespace and comments (threshold 25)\n")


def tidy_output(root: Path, findings: list, canary: bool, extra: str) -> str:
    text = "".join(diagnostic(f"{root}/{path}", 10 + i, name, n) for i, (path, name, n) in enumerate(findings))
    canary_text = diagnostic("/tmp/ctrlpp-function-size.x/all_headers.cpp", 170, "ctrlpp_size_gate_canary", 30)
    return text + (canary_text if canary else "") + extra


def build(root: Path, findings: list, registered: list, canary: bool = True, extra: str = "") -> None:
    git(root, "init", "-q")
    table = "".join(f"| {path} | {name} | {n} | {WHY} |\n" for path, name, n in registered)
    (root / "EXCEPTIONS.md").write_text(FILES + FUNCTIONS + table)
    (root / "tidy.txt").write_text(tidy_output(root.resolve(), findings, canary, extra))
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "fixture")


def run_gate(root: Path, *args: str) -> tuple[int, str]:
    done = subprocess.run([sys.executable, str(GATE), "--source-root", str(root),
                           "--tidy-output", str(root / "tidy.txt"), *args], capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


SOLVE = (HEADER, "solve", 40)
CMAKE = ("lib/probe.cmake", "helper", 30)
CASES = [
    ("no finding at all", [], [], {}, 0, "0 functions over 25 lines; 0 registered"),
    ("an unlisted finding", [SOLVE], [], {}, 1,
     f"error: {HEADER}:solve: 40 lines, over the 25-line ceiling, and is not registered in EXCEPTIONS.md"),
    ("one unlisted beside one registered", [SOLVE, (OTHER, "kernel", 30)], [SOLVE], {}, 1,
     f"error: {OTHER}:kernel: 30 lines"),
    ("a registered function within its allowance", [(HEADER, "solve", 44)], [SOLVE], {}, 0,
     "1 registered functions"),
    ("a registered function past its allowance", [(HEADER, "solve", 45)], [SOLVE], {}, 4,
     f"error: {HEADER}:solve: 45 lines, over the 25-line ceiling, grown from the 40"),
    ("a registered row longer than its function", [(HEADER, "solve", 39)], [SOLVE], {}, 7,
     "records 40 lines but the function is 39"),
    ("a registered row whose function is gone", [], [SOLVE], {}, 3,
     f"error: {HEADER}:solve: the EXCEPTIONS.md row recording 40 lines is stale"),
    ("overloads pair by sorted counts", [(HEADER, "solve", 40), (HEADER, "solve", 30)],
     [(HEADER, "solve", 30), SOLVE], {}, 0, "2 registered functions"),
    ("overloads with one grown", [(HEADER, "solve", 40), (HEADER, "solve", 34)],
     [(HEADER, "solve", 30), SOLVE], {}, 4, "grown from the 30 recorded"),
    ("an overload beyond its rows", [(HEADER, "solve", 40), (HEADER, "solve", 30)], [SOLVE], {}, 1,
     "an overload beyond its rows in EXCEPTIONS.md"),
    ("a CMake function row is the file gate's", [SOLVE], [SOLVE, CMAKE], {}, 0, "1 registered functions"),
    ("no canary in the output", [SOLVE], [SOLVE], {"canary": False}, 6, "canary function"),
    ("a clang-diagnostic-error line", [SOLVE], [SOLVE], {"extra": ERROR}, 6, "could not parse"),
    ("a warning with no line count", [SOLVE], [SOLVE], {"extra": ORPHAN}, 6, "carries no line count"),
    ("a finding outside the source root", [SOLVE], [SOLVE], {"extra": diagnostic("/elsewhere/x.h", 3, "f", 30)},
     6, "outside the checked scope"),
]


def check(label: str, code: int, output: str, want_code: int, want_text: str) -> bool:
    if code == want_code and want_text in output:
        print(f"PASS {label}")
        return True
    print(f"FAIL {label}: exit {code} (wanted {want_code})\n{output.strip()}")
    return False


def case_dir(root: Path, label: str) -> Path:
    case = root / "".join(c if c.isalnum() else "_" for c in label)
    case.mkdir()
    return case


def run_table(root: Path) -> int:
    failures = 0
    for label, findings, registered, options, want_code, want_text in CASES:
        case = case_dir(root, label)
        build(case, findings, registered, **options)
        failures += not check(label, *run_gate(case), want_code, want_text)
    return failures


def run_register_faults(root: Path) -> int:
    case = case_dir(root, "malformed register")
    build(case, [SOLVE], [])
    (case / "EXCEPTIONS.md").write_text(FILES + FUNCTIONS + f"| {HEADER} | solve | a reason |\n")
    failures = not check("a function row with no count", *run_gate(case), 5,
                         "needs a path, a function, a recorded line count and a reason")
    (case / "EXCEPTIONS.md").unlink()
    return failures + (not check("the register itself missing", *run_gate(case), 5, "does not exist"))


def main() -> int:
    with tempfile.TemporaryDirectory() as workspace:
        root = Path(workspace)
        failures = run_table(root) + run_register_faults(root)
    total = len(CASES) + 2
    if failures:
        print(f"{failures} of {total} cases failed")
        return 1
    print(f"{total} cases passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
