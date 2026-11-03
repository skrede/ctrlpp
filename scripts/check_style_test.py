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


NAMING_CASES = [
    ("a private member without m_ fails", lib("x.h", "class probe\n{\n    int count_;\n};\n"),
     14, ("x.h:6: member-naming: probe::count_: a private member takes the m_ prefix; no member name ends in an "
          "underscore",)),
    ("a private member with m_ passes", lib("x.h", "class probe\n{\n    int m_count;\n};\n"), 0, ("no findings",)),
    ("a public aggregate field takes no prefix", lib("x.h", "struct probe\n{\n    int horizon;\n};\n"), 0,
     ("no findings",)),
    ("a protected member without m_ fails", lib("x.h", "class probe\n{\nprotected:\n    int count;\n};\n"), 14,
     ("probe::count: a protected member takes the m_ prefix",)),
    ("a static constexpr member with a trailing underscore fails",
     lib("x.h", "class probe\n{\n    static constexpr int max_phases_ = 11;\n};\n"), 14,
     ("probe::max_phases_: no member name ends in an underscore",)),
    ("a static constexpr member in snake_case passes",
     lib("x.h", "class probe\n{\n    static constexpr int max_phases = 11;\n};\n"), 0, ("no findings",)),
    ("a static inline member with m_ fails", lib("x.h", "class probe\n{\n    static inline int m_count = 0;\n};\n"),
     14, ("probe::m_count: a static member takes no m_ prefix",)),
    ("a public std::function field without _cb fails",
     lib("x.h", "struct probe\n{\n    std::function<void(int, int)> cost;\n};\n"), 14,
     ("probe::cost: a callback member takes the _cb suffix",)),
    ("a public std::function field with _cb passes",
     lib("x.h", "struct probe\n{\n    std::function<void()> cost_cb;\n};\n"), 0, ("no findings",)),
    ("a private optional callback without _cb fails",
     lib("x.h", "class probe\n{\n    std::optional<std::function<void()>> m_hook;\n};\n"), 14,
     ("probe::m_hook: a callback member takes the _cb suffix",)),
    ("a private move-only callback with m_ and _cb passes",
     lib("x.h", "class probe\n{\n    std::move_only_function<void()> m_feed_forward_cb;\n};\n"), 0, ("no findings",)),
    ("a function-pointer member without _cb fails", lib("x.h", "struct probe\n{\n    void (*on_done)(int);\n};\n"),
     14, ("probe::on_done: a callback member takes the _cb suffix",)),
    ("a container of callbacks is not itself a callback",
     lib("x.h", "struct probe\n{\n    std::vector<std::function<void()>> hooks;\n};\n"), 0, ("no findings",)),
    ("each name of a multi-declarator line is reported",
     lib("x.h", "template <typename Scalar>\nclass probe\n{\n    Scalar q0_{}, q1_{};\n    std::array<Scalar, 3> a_, "
                "b_;\n};\n"), 14, ("probe::q0_", "probe::q1_", "probe::a_", "probe::b_")),
    ("function bodies, lambdas and parameters inside a class are not members",
     lib("x.h", "class probe\n{\npublic:\n    void f(int n_ = int{3})\n    {\n        int local_ = n_;\n        auto l = "
                "[](int a_) { int b_ = a_; return b_; };\n        local_ = l(local_);\n    }\n\n    static constexpr "
                "auto make = [](int x_) { return x_; };\n    probe(int c, std::function<void()> f = {}) : m_c{c} {}\n\n"
                "private:\n    int m_c;\n};\n"), 0, ("no findings",)),
    ("an anonymous union's members take the access of their place in the class",
     lib("x.h", "class probe\n{\n    union\n    {\n        int value;\n        float other;\n    };\n};\n"), 14,
     ("probe::<anonymous>::value: a private member takes the m_ prefix",)),
]

INIT_CASES = [
    ("a configuration aggregate keeps its defaults",
     lib("x.h", "struct probe_config\n{\n    int horizon{10};\n    double tol = 1e-9;\n};\n"), 0, ("no findings",)),
    ("the same aggregate named as a result fails",
     lib("x.h", "struct probe_result\n{\n    int horizon{10};\n    double tol = 1e-9;\n};\n"), 13,
     ("x.h:6: member-init: probe_result::horizon is initialized in its declaration", "x.h:7: member-init")),
    ("a configuration struct with a constructor fails",
     lib("x.h", "struct probe_config\n{\n    probe_config() {}\n    int horizon{10};\n};\n"), 13,
     ("probe_config::horizon",)),
    ("a configuration struct with a private member fails",
     lib("x.h", "struct probe_options\n{\n    int horizon{10};\n\nprivate:\n    int m_x;\n};\n"), 13,
     ("probe_options::horizon",)),
    ("a configuration class fails", lib("x.h", "class probe_settings\n{\npublic:\n    int horizon{10};\n};\n"), 13,
     ("probe_settings::horizon",)),
    ("a nested config passes", lib("x.h", "class planner\n{\npublic:\n    struct config\n    {\n        double v_max{1.0};"
                                          "\n    };\n};\n"), 0, ("no findings",)),
    ("a config struct at namespace scope is not the nested carve-out",
     lib("x.h", "struct config\n{\n    double v_max{1.0};\n};\n"), 13, ("config::v_max",)),
    ("a class member initializer fails", lib("x.h", "class probe\n{\n    int m_x{};\n};\n"), 13,
     ("probe::m_x is initialized in its declaration",)),
    ("static constexpr and static inline members are not in-declaration initializers",
     lib("x.h", "class probe\n{\n    static constexpr int n = 3;\n    static inline int k{4};\n};\n"), 0,
     ("no findings",)),
    ("defaulted functions and bit-fields are not members",
     lib("x.h", "struct probe_point\n{\n    probe_point() = default;\n    probe_point& operator=(const probe_point&) = "
                "delete;\n    unsigned ready : 1 {0};\n    double x;\n};\n"), 0, ("no findings",)),
    ("braced default arguments are not member initializers",
     lib("x.h", "class probe\n{\npublic:\n    static int create(int rng = int{}, int likelihood = int{})\n    {\n"
                "        return rng + likelihood;\n    }\n\n    explicit probe(int s = {})\n        : m_s{s}\n    {}\n\n"
                "private:\n    int m_s;\n};\n"), 0, ("no findings",)),
    ("a template argument holding parentheses does not hide an initializer",
     lib("x.h", "template <int N>\nclass probe\n{\n    Eigen::Matrix<double, static_cast<int>(N), 1> m_v{};\n};\n"), 13,
     ("probe::m_v",)),
]

RETURN_CASES = [
    ("a trailing return type fails", lib("x.h", "auto f() -> int;\n"), 12,
     ("x.h:4: return-type: f: a trailing return type counts as auto; spell the return type int in front",)),
    ("a decltype trailing return passes", lib("x.h", "auto f(int x) -> decltype(g(x));\n"), 0, ("no findings",)),
    ("a wrapped head with its arrow on the next line fails",
     lib("x.h", "template <typename T>\ninline auto f(const T& a,\n              const T& b)\n    -> T\n{\n    return a;\n}\n"),
     12, ("x.h:5: return-type: f: a trailing return type",)),
    ("a deduced auto& member fails", lib("x.h", "struct probe\n{\n    constexpr auto& policy();\n};\n"), 12,
     ("probe" if False else "policy: deduced return type; spell it",)),
    ("decltype(auto) fails", lib("x.h", "decltype(auto) get();\n"), 12, ("get: deduced return type",)),
    ("a trailing return on an operator fails",
     lib("x.h", "struct probe\n{\n    auto operator()(int a) const noexcept -> int;\n};\n"), 12,
     ("operator()", "spell the return type int in front")),
    ("an out-of-class definition fails",
     lib("x.h", "template <typename T>\nauto box<T>::get() const -> T\n{\n    return {};\n}\n"), 12, ("box<T>::get",)),
    ("auto variables, auto parameters and lambdas pass",
     lib("x.h", "auto x = 3;\nvoid f(auto x);\ninline auto l = [](int a) -> int { return a; };\nvoid g()\n{\n    auto "
                "h = [](int a) { return a; };\n    auto r = h(1);\n}\n"), 0, ("no findings",)),
]

CALL_CASES = [
    ("a bare fallible call in a body fails", lib("x.h", "void f(int x)\n{\n    create(x);\n}\n"), 19,
     ("x.h:6: discarded-result: the result of create(...) is discarded",)),
    ("a bare call wrapped over two lines fails",
     lib("x.h", "void f(int x)\n{\n    probe::create(x,\n                  x);\n}\n"), 19, ("x.h:6: discarded-result",)),
    ("a qualified bare call fails", lib("x.h", "void f(int a, int b)\n{\n    detail::dare(a, b, a, b);\n}\n"), 19,
     ("dare(...)",)),
    ("a bare call after a condition fails", lib("x.h", "void f(int x)\n{\n    if(x > 0)\n        obj.setup(x);\n}\n"),
     19, ("x.h:7: discarded-result: the result of setup(...)",)),
    ("a leading-return declaration in a class passes",
     lib("x.h", "class pid\n{\npublic:\n    static expected<pid, pid_error> create(const cfg&);\n    static expected<pid, "
                "pid_error>\n        try_create(const cfg&);\n};\n"), 0, ("no findings",)),
    ("a return type on the preceding line at namespace scope passes",
     lib("x.h", "expected<foo, err>\ncreate(int args);\n"), 0, ("no findings",)),
    ("an out-of-class definition head passes",
     lib("x.h", "template <typename T>\nexpected<T, E> type<T>::create(int args)\n{\n    return {};\n}\n"), 0,
     ("no findings",)),
    ("a used result passes",
     lib("x.h", "int f(int x)\n{\n    auto r = create(x);\n    (void)create(x);\n    if(auto s = create(x); s)\n"
                "        return 1;\n    expected<a, b> create(int);\n    update(x);\n    return create(x);\n}\n"), 0,
     ("no findings",)),
    ("a requires-expression is not a body",
     lib("x.h", "template <typename S>\nconcept solver = requires(S s) {\n    s.setup(1);\n};\nvoid f()\n{\n    if "
                "constexpr(requires { place(1); })\n        return;\n}\n"), 0, ("no findings",)),
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


def run_rule_cases(root: Path, rule: str, cases: list) -> int:
    failures = 0
    for index, (label, files, want_code, wanted) in enumerate(cases):
        case = root / f"{rule}-{index}"
        case.mkdir()
        build(case, files)
        code, output = run_lint(case, "--rule", rule)
        failures += not all(check(label, code, output, want_code, text) for text in wanted)
    return failures


def run_list_members(root: Path) -> int:
    case = root / "list-members"
    case.mkdir()
    build(case, lib("x.h", "class outer\n{\n    struct inner\n    {\n        int a{1};\n    };\n    double q0_{}, "
                           "q1_ = 2.0;\n    int m_n;\n};\nstruct outer_config\n{\n    int b{2};\n};\n"))
    code, output = run_lint(case, "--rule", "member-init", "--list-members")
    failures = not all(check("--list-members prints each flagged member with its initializer", code, output, 0, text)
                       for text in (f"{LIB}x.h:8: outer::inner: a {{1}}\n", f"{LIB}x.h:10: outer: q0_ {{}}\n",
                                    f"{LIB}x.h:10: outer: q1_ = 2.0\n"))
    failures += "m_n" in output or "outer_config" in output
    return failures + (not check("--list-members needs --rule member-init", *run_lint(case, "--list-members"), 2,
                                 "--list-members goes with --rule member-init"))


def run_cli(root: Path) -> int:
    case = root / "cli"
    case.mkdir()
    build(case, lib("x.h", "") | {"tests/x.cpp": "\n"})
    failures = not check("a path outside the scope is a usage error", *run_lint(case, "tests/x.cpp"), 2,
                         "outside the lint's scope")
    failures += not check("an unknown rule is a usage error", *run_lint(case, "--rule", "nope"), 2, "invalid choice")
    failures += not check("a given path that cannot be read", *run_lint(case, LIB + "gone.h"), 6, "cannot read")
    failures += not check("--list-rules names every rule and its code", *run_lint(case, "--list-rules"), 0,
                          "include-order 10\nattributes 11\nreturn-type 12\nmember-init 13\nmember-naming 14\n"
                          "pragma-once 15\nguard 16\nclosing-comment 17\nspelling 18\ndiscarded-result 19\n")
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


STRUCTURAL = (("member-naming", NAMING_CASES), ("member-init", INIT_CASES), ("return-type", RETURN_CASES),
              ("discarded-result", CALL_CASES))


def main() -> int:
    groups = (ATTRIBUTE_CASES, TEXT_CASES, INCLUDE_CASES)
    with tempfile.TemporaryDirectory() as workspace:
        root = Path(workspace)
        failures = sum(run_cases(root, cases) for cases in groups[:2])
        failures += run_cases(root, INCLUDE_CASES, {"given paths are linted alone": ("--rule", "include-order",
                                                                                      LIB + "x.h")})
        failures += run_cli(root) + run_in_process()
        for rule, cases in STRUCTURAL:
            failures += run_rule_cases(root, rule, cases)
        failures += run_list_members(root)
    total = sum(len(cases) for cases in groups) + 5 + 4 + sum(len(cases) for _, cases in STRUCTURAL) + 3
    if failures:
        print(f"{failures} of {total} cases failed")
        return 1
    print(f"{total} cases passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
