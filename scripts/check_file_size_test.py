#!/usr/bin/env python3
"""Fixture repositories exercising every register verdict the file-size gate can reach."""

import subprocess
import sys
import tempfile
from pathlib import Path

GATE = Path(__file__).resolve().parent / "check_file_size.py"
WHY = "one cohesive unit whose shared state splitting would scatter across files"
FILES = "# register\n\n## Files\n\n| Path | Lines at registration | Why |\n| --- | --- | --- |\n"
FUNCTIONS = "\n## Functions\n\n| Path | Function | Lines at registration | Why |\n| --- | --- | --- | --- |\n"


def git(root: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(root), "-c", "user.name=fixture", "-c", "user.email=f@x",
                           *args], check=True, capture_output=True, text=True).stdout


def write_tree(root: Path, sources: dict[str, int | str], rows: list, extra: str = "",
               function_rows: str = "") -> None:
    for name, content in sources.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content if isinstance(content, str) else "// x\n" * content)
    table = "".join(f"| {p} | {n} | {why} |\n" for p, n, why in rows)
    (root / "EXCEPTIONS.md").write_text(FILES + table + extra + FUNCTIONS + function_rows)


def build(root: Path, sources: dict[str, int | str], rows: list, extra: str = "",
          function_rows: str = "") -> None:
    git(root, "init", "-q")
    write_tree(root, sources, rows, extra, function_rows)
    git(root, "add", "-A")


def run_gate(root: Path, *args: str) -> tuple[int, str]:
    done = subprocess.run([sys.executable, str(GATE), "--source-root", str(root), *args],
                          capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


def cmake_function(name: str, span: int) -> str:
    return f"function({name} arg)\n" + "  message(STATUS x)\n" * (span - 1) + "endfunction()\n"


BIG = {"lib/big.h": 205}
CASES = [
    ("every overage registered", BIG, [("lib/big.h", 205, WHY)], "", 0, "1 registered files"),
    ("unregistered overage", BIG, [], "", 1, "error: lib/big.h has 205 lines, over the 200-line ceiling, "
     "and is not registered in EXCEPTIONS.md"),
    ("a .cpp overage is seen too", {"lib/big.cpp": 205}, [], "", 1, "lib/big.cpp"),
    ("C++ files outside lib are out of scope",
     {"tests/big.cpp": 900, "benchmarks/big.h": 900, "src/big.h": 900, "lib/small.h": 10}, [], "", 0,
     "checked 1 tracked files (0 CMake)"),
    ("boundary: the ceiling itself passes", {"lib/edge.h": 200}, [], "", 0, "0 registered"),
    ("boundary: one line over fails", {"lib/edge.h": 201}, [], "", 1, "has 201 lines"),
    ("stale row, file under the ceiling",
     {"lib/small.h": 10, "lib/big.h": 205}, [("lib/small.h", 205, WHY), ("lib/big.h", 205, WHY)],
     "", 3, "error: lib/small.h: the EXCEPTIONS.md row is stale"),
    ("stale row, no such file", BIG, [("lib/big.h", 205, WHY), ("lib/gone.h", 205, WHY)], "",
     3, "error: lib/gone.h: the EXCEPTIONS.md row is stale"),
    ("grown past the allowance", {"lib/big.h": 227}, [("lib/big.h", 205, WHY)], "", 4,
     "error: lib/big.h has grown to"),
    ("at exactly the allowance passes", {"lib/big.h": 226}, [("lib/big.h", 205, WHY)], "", 0, "1 reg"),
    ("recorded count above the file's length", BIG, [("lib/big.h", 300, WHY)], "",
     7, "records 300 lines but the file has 205"),
    ("a row with no count column", BIG, [], "| lib/big.h | a reason |\n", 5,
     "needs a path, a recorded line count and a reason"),
    ("a row with an empty reason", BIG, [("lib/big.h", 205, "")], "", 5, "reason is empty"),
    ("a reason restating the count", BIG, [("lib/big.h", 205, "it is 205 lines")],
     "", 5, "restates the unit's line count"),
    ("a reason naming no split harm", BIG, [("lib/big.h", 205, "a large file")],
     "", 5, "must name what splitting"),
    ("a reason promising a future split", BIG,
     [("lib/big.h", 205, "cohesive today; splitting is planned")], "", 5, "promises future work"),
    ("a path outside the scope", BIG,
     [("lib/big.h", 205, WHY), ("tests/big.cpp", 205, WHY)], "", 5, "outside the checked scope"),
    ("a path registered twice", BIG, [("lib/big.h", 205, WHY), ("lib/big.h", 205, WHY)], "", 5,
     "registered twice"),
    ("a row under an unknown heading", BIG, [("lib/big.h", 205, WHY)],
     f"\n## Notes\n\n| lib/other.h | 205 | {WHY} |\n", 5, "neither 'Files' nor 'Functions'"),
    ("a path in backticks", BIG, [("`lib/big.h`", 205, WHY)], "", 0, "1 registered"),
    ("a path carrying a space", {"lib/my file.h": 205}, [("lib/my file.h", 205, WHY)], "", 0, "1 reg"),
    ("no file in scope anywhere", {"README.md": 5, "src/a.h": 300}, [], "", 6, "no tracked file in"),
    ("an over-limit CMake file under lib", {"lib/probe.cmake": 201}, [], "", 1,
     "error: lib/probe.cmake has 201 lines"),
    ("an over-limit root CMakeLists.txt", {"CMakeLists.txt": 201}, [], "", 1,
     "error: CMakeLists.txt has 201 lines"),
    ("an over-limit .cmake.in", {"cmake/config.cmake.in": 201}, [], "", 1, "cmake/config.cmake.in"),
    ("an over-limit tests/CMakeLists.txt is a guideline",
     {"tests/CMakeLists.txt": 201, "examples/x.cmake": 201, "benchmarks/CMakeLists.txt": 201,
      "lib/small.h": 10}, [], "", 0, "checked 1 tracked files (0 CMake)"),
    ("a registered CMake file", {"lib/probe.cmake": 201}, [("lib/probe.cmake", 201, WHY)], "", 0,
     "checked 1 tracked files (1 CMake); 1 registered files"),
    ("a CMake function at the ceiling", {"lib/probe.cmake": cmake_function("ok", 25)}, [], "", 0,
     "0 registered functions"),
    ("a CMake function one line over", {"lib/probe.cmake": cmake_function("long_one", 26)}, [], "", 1,
     "error: lib/probe.cmake:long_one: 26 lines, over the 25-line ceiling, and is not registered"),
    ("a CMake macro one line over", {"lib/probe.cmake": cmake_function("x", 26).replace("function", "macro")},
     [], "", 1, "error: lib/probe.cmake:x: 26 lines"),
    ("an upper-case CMake function", {"CMakeLists.txt": cmake_function("UP", 26).upper()}, [], "", 1,
     "error: CMakeLists.txt:UP: 26 lines"),
    ("a long CMake function in tests is a guideline",
     {"tests/CMakeLists.txt": cmake_function("long_one", 40), "lib/small.h": 10}, [], "", 0, "checked 1"),
]

CMAKE_ROW = f"| lib/probe.cmake | long_one | 26 | {WHY} |\n"
CPP_ROW = f"| lib/big.h | solve | 40 | {WHY} |\n"
FUNCTION_CASES = [
    ("a registered CMake function", 26, CMAKE_ROW, 0, "1 registered functions"),
    ("a registered CMake function within its allowance", 29, CMAKE_ROW, 0, "1 registered functions"),
    ("a registered CMake function past its allowance", 30, CMAKE_ROW, 4, "past its allowance of 29"),
    ("a registered CMake function shrunk under its row", 25, CMAKE_ROW, 3,
     "error: lib/probe.cmake:long_one: the EXCEPTIONS.md row recording 26 lines is stale"),
    ("a CMake row longer than its function", 27, CMAKE_ROW.replace("| 26 |", "| 28 |"), 7,
     "records 28 lines but the function is 27"),
    ("a C++ function row is the function gate's", 26, CMAKE_ROW + CPP_ROW, 0, "1 registered functions"),
    ("a function row with no count", 26, "| lib/probe.cmake | long_one | a reason |\n", 5,
     "needs a path, a function, a recorded line count and a reason"),
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
    for label, sources, rows, extra, want_code, want_text in CASES:
        case = case_dir(root, label)
        build(case, sources, rows, extra)
        failures += not check(label, *run_gate(case), want_code, want_text)
    return failures


def run_function_rows(root: Path) -> int:
    failures = 0
    for label, span, function_rows, want_code, want_text in FUNCTION_CASES:
        case = case_dir(root, label)
        build(case, {"lib/probe.cmake": cmake_function("long_one", span)}, [], "", function_rows)
        failures += not check(label, *run_gate(case), want_code, want_text)
    return failures


def run_tree_faults(root: Path) -> int:
    absent = case_dir(root, "absent")
    build(absent, {"lib/big.h": 205, "lib/kept.h": 10}, [("lib/big.h", 205, WHY)])
    (absent / "lib/big.h").unlink()
    failures = not check("a tracked file missing from the working tree", *run_gate(absent),
                         6, "error: lib/big.h is tracked but absent from the working tree")
    bare = case_dir(root, "no register")
    build(bare, BIG, [])
    (bare / "EXCEPTIONS.md").unlink()
    return failures + (not check("the register itself missing", *run_gate(bare), 5, "does not exist"))


def run_counting(root: Path) -> int:
    tail = case_dir(root, "no trailing newline")
    build(tail, {"lib/big.h": 200}, [])
    (tail / "lib/big.h").write_text("// x\n" * 200 + "// last")
    return not check("a final line with no newline still counts", *run_gate(tail), 1, "has 201 lines")


def main() -> int:
    with tempfile.TemporaryDirectory() as workspace:
        root = Path(workspace)
        failures = run_table(root) + run_function_rows(root) + run_tree_faults(root) + run_counting(root)
    total = len(CASES) + len(FUNCTION_CASES) + 3
    if failures:
        print(f"{failures} of {total} cases failed")
        return 1
    print(f"{total} cases passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
