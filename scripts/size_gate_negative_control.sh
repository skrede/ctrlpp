#!/usr/bin/env bash
# Prove that the size gates are able to fail: on an over-limit file, CMake file, CMake function or C++
# function nobody registered, and on every way a register row can stop describing its unit.
#
# The gates measure tracked files only, so an untracked plant would be invisible to them. Every plant
# therefore goes into a throwaway repository holding a copy of the real tree, and the real working
# tree is never written. The copy registers the overages it starts with, so it starts green and each
# plant is the only thing that can turn it red.
#
# The scope names the gates to prove: files, functions, or all (the default). The function plants run a
# real clang-tidy, found through CLANG_TIDY, then clang-tidy-18, then the llvm18 install path. Without
# one, scope all skips them visibly and scope functions fails. CTRLPP_TIDY_EXTRA_ARG reaches the
# function gate unchanged; it must carry the include path of the argmin headers, which the library
# includes unconditionally, and where needed --gcc-install-dir=... for a libstdc++ clang-18 cannot parse.
#
#   scripts/size_gate_negative_control.sh [files|functions|all]
set -euo pipefail

SCOPE="${1:-all}"
case "$SCOPE" in
    files|functions|all) ;;
    *) echo "usage: $0 [files|functions|all]" >&2; exit 2 ;;
esac

cd "$(dirname "$0")/.."

PROBE="lib/ctrlpp/include/ctrlpp/size_gate_probe.h"
REGISTERED_PROBE="lib/ctrlpp/include/ctrlpp/size_gate_registered.h"
CMAKE_PROBE="lib/probe.cmake"
GUIDELINE_PROBE="tests/probe/CMakeLists.txt"
REGISTER="EXCEPTIONS.md"
WORK=""

cleanup() { if [ -n "$WORK" ]; then rm -rf "$WORK"; fi; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

failures=0
pass() { printf '  PASS  %s\n' "$1"; }
fail() { printf '  FAIL  %s\n' "$1"; failures=$((failures + 1)); }

in_copy() { git -C "$WORK" -c user.name=control -c user.email=control@localhost \
    -c commit.gpgsign=false -c core.hooksPath=/dev/null "$@"; }

gate_status() {
    local rc=0
    (cd "$WORK" && "$@") >/dev/null 2>&1 || rc=$?
    echo "$rc"
}

expect() {
    local want="$1" label="$2" got
    shift 2
    got="$(gate_status "$@")"
    if [ "$got" -eq "$want" ]; then pass "$label (exit $got)"; else fail "$label: exit $got, wanted $want"; fi
}

plant_lines() {
    mkdir -p "$(dirname "$WORK/$1")"
    for _ in $(seq "$2"); do echo "$3"; done > "$WORK/$1"
    in_copy add "$1"
}

unplant() { in_copy rm -q --cached "$1"; rm -f "$WORK/$1"; }
restore() { in_copy checkout -q -- "$@"; }

plant_function() {
    {
        echo '#ifndef HPP_GUARD_CTRLPP_SIZE_GATE_PROBE_H'
        echo '#define HPP_GUARD_CTRLPP_SIZE_GATE_PROBE_H'
        echo 'inline void size_gate_probe()'
        echo '{'
        for _ in $(seq 29); do echo '    static_cast<void>(0);'; done
        echo '}'
        echo '#endif'
    } > "$WORK/$PROBE"
    in_copy add "$PROBE"
}

plant_cmake_function() {
    {
        echo 'function(size_gate_probe)'
        for _ in $(seq 25); do echo '  message(STATUS probe)'; done
        echo 'endfunction()'
    } > "$WORK/$CMAKE_PROBE"
    in_copy add "$CMAKE_PROBE"
}

find_tidy() {
    local candidate
    for candidate in "${CLANG_TIDY:-}" clang-tidy-18 /usr/lib/llvm18/bin/clang-tidy; do
        if [ -n "$candidate" ] && command -v "$candidate" >/dev/null 2>&1; then
            command -v "$candidate"
            return
        fi
    done
}

edit_register() {
    python3 - "$WORK/$REGISTER" "$@" <<'PY'
import re
import sys

path, action, *rest = sys.argv[1:]
lines = open(path, encoding="utf-8").read().splitlines()
WHY = "a throwaway registration in the control's copy, whose splitting is not under test"
FILE_FINDING = re.compile(r"^error: (\S.*?) has (\d+) lines, over the 200-line ceiling, and is not registered")
FUNCTION_FINDING = re.compile(r"^error: ([^:]+):(.*): (\d+) lines, over the 25-line ceiling, and is not registered")


def insert(section, rows):
    start = lines.index(f"## {section}")
    end = start + 1
    while end < len(lines) and not lines[end].startswith("#"):
        end += 1
    last = max(i for i in range(start, end) if lines[i].startswith("|"))
    lines[last + 1:last + 1] = rows


def rows_of(section):
    start = lines.index(f"## {section}")
    return [i for i in range(start + 1, len(lines)) if lines[i].startswith("| lib/")
            and not any(line.startswith("#") for line in lines[start + 1:i])]


if action == "register":
    findings = open(rest[0], encoding="utf-8").read().splitlines()
    files = [m for m in map(FILE_FINDING.match, findings) if m]
    functions = [m for m in map(FUNCTION_FINDING.match, findings) if m]
    insert("Files", [f"| {m[1]} | {m[2]} | {WHY} |" for m in files])
    insert("Functions", [f"| {m[1]} | {m[2]} | {m[3]} | {WHY} |" for m in functions])
elif action == "row":
    insert(rest[0], [rest[1]])
elif action == "first-file":
    print(lines[rows_of("Files")[0]].split("|")[1].strip(), lines[rows_of("Files")[0]].split("|")[2].strip())
    sys.exit(0)
elif action == "inflate-file":
    index = rows_of("Files")[0]
    cells = lines[index].split("|")
    cells[2] = f" {int(cells[2]) + 1} "
    lines[index] = "|".join(cells)
elif action == "lower-function":
    if not rows_of("Functions"):
        sys.exit(1)
    index = max(rows_of("Functions"), key=lambda i: int(lines[i].split("|")[3]))
    cells = lines[index].split("|")
    cells[3] = f" {int(cells[3]) * 10 // 11 - 1} "
    lines[index] = "|".join(cells)
open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")
PY
}

before="$(git status --porcelain)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ctrlpp-size-gate.XXXXXX")"
git ls-files -z -- lib scripts EXCEPTIONS.md CMakeLists.txt cmake validation \
    | tar --null -T - -cf - | tar -xf - -C "$WORK"
in_copy init -q
in_copy add -A
# The real register may hold no file row, so the copy always carries one over-limit file to register.
plant_lines "$REGISTERED_PROBE" 220 '// probe'

FILE_GATE=(python3 scripts/check_file_size.py)
TIDY="$(find_tidy)"
FUNCTION_GATE=(python3 scripts/check_function_size.py --clang-tidy "$TIDY")
RUN_FUNCTIONS=0
if [ "$SCOPE" != files ] && [ -n "$TIDY" ]; then RUN_FUNCTIONS=1; fi

collect_overages() {
    local rc=0
    (cd "$WORK" && "$@") 2> "$WORK/.gate" >/dev/null || rc=$?
    if [ "$rc" -gt 1 ]; then
        fail "the copy cannot be measured, so its overages cannot be registered: $* exited $rc"
        tail -n 5 "$WORK/.gate" | sed 's/^/        /'
    fi
    cat "$WORK/.gate" >> "$WORK/.overages"
    rm -f "$WORK/.gate"
}

: > "$WORK/.overages"
collect_overages "${FILE_GATE[@]}"
if [ "$RUN_FUNCTIONS" -eq 1 ]; then collect_overages "${FUNCTION_GATE[@]}"; fi
edit_register register "$WORK/.overages"
rm -f "$WORK/.overages"
in_copy add "$REGISTER"
in_copy commit -q -m "copy of the real tree with its overages registered"

prove_file_gate() {
    echo "proving the file-size gate can fail:"
    local recorded count
    read -r recorded count <<< "$(edit_register first-file)"

    expect 0 "the untouched copy passes" "${FILE_GATE[@]}"

    plant_lines "$PROBE" 201 '// probe'
    expect 1 "a planted 201-line lib/ header in no register is caught" "${FILE_GATE[@]}"
    unplant "$PROBE"

    for _ in $(seq $(( count / 10 + 2 ))); do echo '// probe'; done >> "$WORK/$recorded"
    expect 4 "a registered file grown past its allowance is caught" "${FILE_GATE[@]}"
    restore "$recorded"

    edit_register row Files "| $PROBE | 201 | a throwaway row whose splitting is not under test |"
    expect 3 "a stale row for a file that does not exist is caught" "${FILE_GATE[@]}"
    restore "$REGISTER"

    edit_register inflate-file
    expect 7 "a row recording more lines than its file has is caught" "${FILE_GATE[@]}"
    restore "$REGISTER"

    edit_register row Files "| $recorded | a reason |"
    expect 5 "a malformed row is caught" "${FILE_GATE[@]}"
    restore "$REGISTER"

    plant_lines "$CMAKE_PROBE" 201 '# probe'
    expect 1 "a planted 201-line lib/ CMake file in no register is caught" "${FILE_GATE[@]}"
    unplant "$CMAKE_PROBE"

    plant_cmake_function
    expect 1 "a planted 26-line CMake function in no register is caught" "${FILE_GATE[@]}"
    unplant "$CMAKE_PROBE"

    plant_lines "$GUIDELINE_PROBE" 201 '# probe'
    expect 0 "a 201-line tests/ CMake file stays a guideline" "${FILE_GATE[@]}"
    unplant "$GUIDELINE_PROBE"
}

prove_function_gate() {
    echo "proving the function-size gate can fail:"
    if [ -z "$TIDY" ] && [ "$SCOPE" = functions ]; then
        fail "function plants: no clang-tidy found; set CLANG_TIDY"
        return
    fi
    if [ -z "$TIDY" ]; then
        echo "  SKIP  function plants: no clang-tidy found; set CLANG_TIDY to run them"
        return
    fi
    expect 0 "the untouched copy passes the function gate" "${FUNCTION_GATE[@]}"

    plant_function
    expect 1 "a planted 30-line function in no register is caught by $TIDY" "${FUNCTION_GATE[@]}"
    unplant "$PROBE"

    if edit_register lower-function; then
        expect 4 "a registered function past its lowered row's allowance is caught" "${FUNCTION_GATE[@]}"
        restore "$REGISTER"
    else
        fail "no registered function to lower; the copy measured no function over the ceiling"
    fi
}

if [ "$SCOPE" != functions ]; then prove_file_gate; fi
if [ "$SCOPE" != files ]; then prove_function_gate; fi

if [ "$(git status --porcelain)" = "$before" ]; then
    pass "the real working tree is unchanged"
else
    fail "the real working tree changed while the control ran"
fi

echo
if [ "$failures" -ne 0 ]; then
    echo "the size gates are NOT proven able to fail: $failures of the control's properties do not hold"
    exit 1
fi

echo "the size gates discriminate: the copy with its overages registered passes, and every planted"
echo "overage and every row that no longer describes its unit goes red"
