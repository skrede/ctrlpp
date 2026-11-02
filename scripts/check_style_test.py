#!/usr/bin/env python3
"""Fixture repositories proving each style rule fails on a minimal violation and passes on its known traps."""

import re
import subprocess
import sys
import tempfile
from pathlib import Path

import style_rules_text
from style_source import Source, blank

LINT = Path(__file__).resolve().parent / "check_style.py"
LIB = "lib/ctrlpp/include/ctrlpp/"
ESP = "examples/embedded/esp32/main/"


def guard_of(path: str) -> str:
    prefix, under = ("HPP_GUARD_", path[len("lib/ctrlpp/include/"):]) if path.startswith("lib/") else \
        ("HPP_GUARD_CTRLPP_EXAMPLES_", path[len("examples/"):])
    return prefix + re.sub(r"[^A-Za-z0-9]", "_", under).upper()


def header(path: str, body: str) -> str:
    guard = guard_of(path)
    return f"#ifndef {guard}\n#define {guard}\n\n{body}\n#endif\n"


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


def build(root: Path, files: dict[str, str]) -> None:
    git(root, "init", "-q")
    for name, content in files.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    if files:
        git(root, "add", "-A")


def run_lint(root: Path, *args: str) -> tuple[int, str]:
    done = subprocess.run([sys.executable, str(LINT), "--source-root", str(root), *args],
                          capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


def lib(name: str, body: str) -> dict[str, str]:
    return {LIB + name: header(LIB + name, body)}


ATTRIBUTE_CASES = [
    ("a class-level [[nodiscard]] fails", lib("x.h", "class [[nodiscard]] x {};\n"), 11,
     "x.h:4: attributes: [[nodiscard]]"),
    ("the attribute inside a line comment and a string passes",
     lib("x.h", '// carries [[nodiscard]]\nconst char* s = "[[nodiscard]]";\n'), 0, "no findings"),
    ("the attribute inside a block comment and a raw string passes",
     lib("x.h", '/* [[nodiscard]]\n */\nconst char* s = R"x( )" [[nodiscard]] )x";\n'), 0, "no findings"),
    ("a digit separator does not open a literal that hides an attribute",
     lib("x.h", "int n = 1'000'000;\n[[nodiscard]] int f();\n"), 11, "x.h:5: attributes"),
    ("a character literal of a bracket does not hide an attribute",
     lib("x.h", "char c = '[';\n[[deprecated]] int f();\n"), 11, "x.h:5: attributes: [[deprecated]]"),
    ("a lambda used as a subscript passes", lib("x.h", "int g(int* a) { return a[[&]{ return 0; }()]; }\n"), 0,
     "no findings"),
    ("the three load-bearing attributes pass",
     lib("x.h", "[[noreturn]] void f();\n[[maybe_unused]] int x;\nvoid g(int k) { switch(k) { case 0: "
                "[[fallthrough]]; default: break; } }\n"), 0, "no findings"),
    ("one disallowed attribute in a list fails", lib("x.h", "[[maybe_unused, nodiscard]] int f();\n"), 11,
     "[[nodiscard]]"),
    ("gnu::noinline is carved out in the embedded examples", {ESP + "probe.cpp": "[[gnu::noinline]] void f() {}\n"},
     0, "no findings"),
    ("gnu::noinline in the library fails", lib("x.h", "[[gnu::noinline]] void f() {}\n"), 11, "[[gnu::noinline]]"),
    ("files outside lib/ and examples/ are not linted", {"tests/x.cpp": "[[nodiscard]] int f();\n",
                                                          LIB + "ok.h": header(LIB + "ok.h", "")}, 0,
     "checked 1 files, no findings"),
]

TEXT_CASES = [
    ("#pragma once fails", {LIB + "x.h": "#pragma once\n"}, 15, "x.h:1: pragma-once"),
    ("#pragma once inside a comment passes", lib("x.h", "// #pragma once\n"), 0, "no findings"),
    ("a misnamed guard fails", {LIB + "control/x.h": "#ifndef X_H\n#define X_H\n#endif\n"}, 16,
     "the guard is X_H; it should be HPP_GUARD_CTRLPP_CONTROL_X_H"),
    ("a header with no guard fails", {LIB + "x.h": "int f();\n"}, 16, "x.h:1: guard"),
    ("code after the guard's #endif fails", {LIB + "x.h": header(LIB + "x.h", "") + "int f();\n"}, 16, "guard"),
    ("a .cpp file without a guard passes", {"examples/x.cpp": "int main() {}\n"}, 0, "no findings"),
    ("a correct library guard passes", {LIB + "control/pid.h": "// a\n#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H\n"
                                        "#define HPP_GUARD_CTRLPP_CONTROL_PID_H\n#endif\n"}, 0, "no findings"),
    ("a correct examples guard passes",
     {ESP + "cycle_record.h": "#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_CYCLE_RECORD_H\n"
                              "#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_CYCLE_RECORD_H\n#endif\n"}, 0,
     "no findings"),
    ("a closing namespace comment fails", lib("x.h", "namespace ctrlpp\n{\n} // namespace ctrlpp\n"), 17,
     "x.h:6: closing-comment"),
    ("a block-comment namespace closer fails", {"examples/x.cpp": "namespace {\n} /* namespace */\n"}, 17,
     "x.cpp:2: closing-comment"),
    ("a comment after the guard's #endif fails",
     {LIB + "x.h": header(LIB + "x.h", "").replace("#endif\n", "#endif // HPP_GUARD_CTRLPP_X_H\n")}, 17,
     "x.h:5: closing-comment: a comment after the guard's #endif"),
    ("a conditional's #endif comment passes",
     lib("x.h", "#if defined(A)\nint f();\n#endif // A && __has_include(<b.h>)\n"), 0, "no findings"),
    ("a namespace closer quoted inside a comment passes", {"examples/x.cpp": "/* } // namespace */\n"}, 0,
     "no findings"),
    ("a British word fails", lib("x.h", "// the behaviour when optimising\n"), 18,
     "'behaviour' is British spelling; write 'behavior'"),
    ("the -ise forms name the American form", {"examples/x.cpp": "// Centralises the normalisation\n"}, 18,
     "write 'Centralizes'"),
    ("British spelling in a CMake file fails", {"cmake/x.cmake": "# colour the output\n",
                                                "CMakeLists.txt": "# fine\n"}, 18, "x.cmake:1: spelling"),
    ("analogue is a US variant and passes", {"examples/x.cpp": "// the analogue of the discrete filter\n"}, 0,
     "no findings"),
    ("a British word inside a URL passes",
     {"examples/x.cpp": "// https://example.org/optimisation-behaviour.html\n"}, 0, "no findings"),
    ("American spellings pass", {"examples/x.cpp": "// analyzes the modeled behavior; a precise exercise\n"}, 0,
     "no findings"),
]


def includes(*lines: str) -> str:
    return "\n".join(lines) + "\n"


PID_GAINS = includes('#include "ctrlpp/control/detail/pid_config_checks.h"', "",
                     '#include "ctrlpp/control/pid_error.h"', '#include "ctrlpp/control/pid_config.h"',
                     '#include "ctrlpp/control/pid_policies.h"', "", '#include "ctrlpp/types.h"',
                     '#include "ctrlpp/expected.h"', "", "#include <cmath>", "#include <limits>",
                     "#include <cstddef>", "", "namespace ctrlpp::detail\n{\n}")
TREE = {LIB + name: header(LIB + name, "") for name in (
    "types.h", "expected.h", "control/pid_error.h", "control/pid_config.h", "control/pid_policies.h",
    "control/detail/pid_config_checks.h", "sysid/rls.h", "util/concepts.h")}
ESP_TREE = {path: header(path, "") for path in (ESP + "eigen_alloc_sentinel.h", ESP + "posture.h", ESP + "telemetry.h",
                                                "examples/embedded/shared/dsp_demo.h")}


def tree(path: str, block: str) -> dict[str, str]:
    files = TREE | ESP_TREE
    files[path] = header(path, block) if path.endswith(".h") else block
    return files


INCLUDE_CASES = [
    ("a detail file with its own detail folder first passes", tree(LIB + "control/detail/pid_gains.h", PID_GAINS), 0,
     "no findings"),
    ("the root folder ahead of the own folder fails",
     tree(LIB + "sysid/recursive_arx.h", includes('#include "ctrlpp/types.h"', '#include "ctrlpp/expected.h"', "",
                                                  '#include "ctrlpp/sysid/rls.h"', "", '#include "ctrlpp/util/concepts.h"')),
     10, "recursive_arx.h:7: include-order: the file's own folder ctrlpp/sysid must lead"),
    ("a detail file whose module folder precedes its detail folder fails",
     tree(LIB + "control/detail/pid_gains.h", includes('#include "ctrlpp/control/pid_error.h"', "",
                                                       '#include "ctrlpp/control/detail/pid_config_checks.h"')),
     10, "own folder ctrlpp/control/detail must lead the project includes, ahead of ctrlpp/control"),
    ("ESP-IDF quoted headers are third-party and the required orders pass",
     tree(ESP + "app_main.cpp", includes("// Armed first.", '#include "eigen_alloc_sentinel.h"', "", "// About it.",
                                         "//", "// More.", "", '#include "posture.h"', '#include "telemetry.h"', "",
                                         '#include "dsp_demo.h"', "", '#include "freertos/FreeRTOS.h"',
                                         '#include "freertos/task.h"', "", '#include "esp_log.h"', "",
                                         "#include <cstdint>", "", "namespace {}")), 0, "no findings"),
    ("the FreeRTOS order outside its listed path still fails",
     tree("examples/other/main.cpp", includes('#include "freertos/FreeRTOS.h"', '#include "freertos/task.h"')),
     10, '"freertos/task.h" sorts before "freertos/FreeRTOS.h"'),
    ("an unlisted order beside a listed one still fails",
     tree(ESP + "app_main.cpp", includes('#include "freertos/FreeRTOS.h"', '#include "freertos/task.h"', "",
                                         '#include "esp_system.h"', '#include "esp_log.h"')),
     10, '"esp_log.h" sorts before "esp_system.h"'),
    ("a listed include sharing a group with another folder still fails",
     tree(ESP + "app_main.cpp", includes('#include "eigen_alloc_sentinel.h"', '#include "dsp_demo.h"')),
     10, '"dsp_demo.h" (examples/embedded/shared) shares a group with "eigen_alloc_sentinel.h"'),
    ("a conditional region is checked on its own",
     tree(LIB + "x.h", includes('#include "ctrlpp/types.h"', "", "#include <memory>", "#include <version>", "",
                                "#if defined(__cpp_lib_expected)", "    #include <expected>", "#endif", "",
                                "#if defined(__GLIBCXX__)", "#include <cxxabi.h>", "#endif", "", "namespace x {}")),
     0, "no findings"),
    ("a forwarding header passes",
     {LIB + "pid.h": "// ctrlpp/pid.h\n#ifndef HPP_GUARD_CTRLPP_PID_H\n#define HPP_GUARD_CTRLPP_PID_H\n"
                     '#include "ctrlpp/control/pid.h"\n#endif\n', LIB + "control/pid.h": header(LIB + "control/pid.h", "")},
     0, "no findings"),
    ("two blank lines between groups fail",
     tree(LIB + "x.h", includes('#include "ctrlpp/types.h"', "", "", "#include <cmath>")), 10,
     "2 blank lines before this line"),
    ("two blank lines before a conditional fail",
     tree(LIB + "x.h", includes("#include <cmath>", "", "", "#if A", "#include <cstdio>", "#endif")), 10,
     "2 blank lines"),
    ("a standard header ahead of a third-party one fails",
     tree(LIB + "x.h", includes("#include <cmath>", "", "#include <Eigen/Dense>")), 10,
     "<Eigen/Dense> is a third-party include after the standard group of <cmath>"),
    ("two folders in one group fail",
     tree(LIB + "x.h", includes("#include <Eigen/Dense>", "#include <unsupported/Eigen/MatrixFunctions>")), 10,
     "shares a group with <Eigen/Dense> (third-party Eigen)"),
    ("a folder split across groups fails",
     tree(LIB + "x.h", includes("#include <cmath>", "", "#include <cstdio>")), 10,
     "the standard library is split into more than one group"),
    ("a longer name ahead of a shorter one fails",
     tree(LIB + "x.h", includes("#include <cstddef>", "#include <cmath>")), 10, "<cmath> sorts before <cstddef>"),
    ("equal lengths sort alphabetically", tree(LIB + "x.h", includes("#include <cmath>", "#include <array>")), 10,
     "<array> sorts before <cmath>"),
    ("a commented-out include is not part of the block",
     tree(LIB + "x.h", includes("#include <cmath>", "// #include <a>", "/* #include <b> */", "#include <cstdio>")),
     0, "no findings"),
    ("given paths are linted alone", tree(LIB + "x.h", includes("#include <cstddef>", "#include <cmath>")) |
     {LIB + "y.h": header(LIB + "y.h", "#include <cstddef>\n#include <cmath>\n")}, 10, "checked 1 files"),
]


def check(label: str, code: int, output: str, want_code: int, want_text: str) -> bool:
    if code == want_code and want_text in output:
        return True
    print(f"FAIL {label}: exit {code} (want {want_code}), expected {want_text!r} in:\n{output}")
    return False


def run_cases(root: Path, cases: list, extra: dict | None = None) -> int:
    failures = 0
    for index, (label, files, want_code, want_text) in enumerate(cases):
        case = root / f"{id(cases)}-{index}"
        case.mkdir()
        build(case, files)
        args = (extra or {}).get(label, ())
        failures += not check(label, *run_lint(case, *args), want_code, want_text)
    return failures


def run_cli(root: Path) -> int:
    case = root / "cli"
    case.mkdir()
    build(case, lib("x.h", "") | {"tests/x.cpp": "\n"})
    failures = not check("a path outside the scope is a usage error", *run_lint(case, "tests/x.cpp"), 2,
                         "outside the lint's scope")
    failures += not check("an unknown rule is a usage error", *run_lint(case, "--rule", "nope"), 2, "invalid choice")
    failures += not check("a given path that cannot be read", *run_lint(case, LIB + "gone.h"), 6, "cannot read")
    failures += not check("--list-rules names every rule and its code", *run_lint(case, "--list-rules"), 0,
                          "include-order 10\nattributes 11\npragma-once 15\nguard 16\nclosing-comment 17\nspelling 18")
    empty = root / "empty"
    empty.mkdir()
    build(empty, {})
    return failures + (not check("an empty repository has nothing to check", *run_lint(empty), 6, "nothing to check"))


def run_in_process() -> int:
    text = 'a = 1\'0; b = u8"x\\"y"; c = R"q(")q"; /* e\n f */ g // h\n'
    blanked = blank(text)
    failures = not check("blanking keeps literal delimiters", 0, blanked, 0, 'a = 1\'0; b = u8"    "; c = R"     "; ')
    failures += not check("blanking keeps newline positions", 0, str([i for i, c in enumerate(blanked) if c == "\n"]),
                          0, str([i for i, c in enumerate(text) if c == "\n"]))
    failures += not check("blanking spans a block comment's newline", 0, blanked, 0, '\n      g     \n')
    style_rules_text.ALLOWED_PHRASES = (("Optimisation in Practice", "a cited title"),)
    found = style_rules_text.check_spelling(Source("examples/x.cpp", "// see Optimisation in Practice\n", "", frozenset()))
    style_rules_text.ALLOWED_PHRASES = ()
    return failures + (not check("an allowed quoted title passes", 0, str(found), 0, "[]"))


def main() -> int:
    groups = (ATTRIBUTE_CASES, TEXT_CASES, INCLUDE_CASES)
    with tempfile.TemporaryDirectory() as workspace:
        root = Path(workspace)
        failures = sum(run_cases(root, cases) for cases in groups[:2])
        failures += run_cases(root, INCLUDE_CASES, {"given paths are linted alone": ("--rule", "include-order",
                                                                                      LIB + "x.h")})
        failures += run_cli(root) + run_in_process()
    total = sum(len(cases) for cases in groups) + 5 + 4
    if failures:
        print(f"{failures} of {total} cases failed")
        return 1
    print(f"{total} cases passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
