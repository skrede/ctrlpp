#!/usr/bin/env bash
# Prove that every style-lint rule is able to fail: a probe header that passes all eleven rules is planted
# into a throwaway copy of the tree, then rewritten with one violation per rule, and each rewrite must
# fail that rule with its own exit code.
#
# The lint reads tracked files only, so the probe goes into a fresh repository holding a copy of the
# real tree, and the real working tree is never written.
#
#   scripts/style_lint_negative_control.sh
set -euo pipefail

cd "$(dirname "$0")/.."

PROBE="lib/ctrlpp/include/ctrlpp/style_lint_probe.h"
WORK=""

cleanup() { if [ -n "$WORK" ]; then rm -rf "$WORK"; fi; }
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

failures=0
pass() { printf '  PASS  %s\n' "$1"; }
fail() { printf '  FAIL  %s\n' "$1"; failures=$((failures + 1)); }

in_copy() { git -C "$WORK" -c core.hooksPath=/dev/null "$@"; }

lint_status() {
    local rc=0
    (cd "$WORK" && python3 scripts/check_style.py "$@" "$PROBE") >/dev/null 2>&1 || rc=$?
    echo "$rc"
}

write_probe() {
    cat > "$WORK/$PROBE" <<'CPP'
#ifndef HPP_GUARD_CTRLPP_STYLE_LINT_PROBE_H
#define HPP_GUARD_CTRLPP_STYLE_LINT_PROBE_H

#include "ctrlpp/types.h"

#include <array>
#include <cstdint>

namespace ctrlpp
{

struct probe_config
{
    std::int32_t horizon{10};
};

class style_lint_probe
{
public:
    explicit style_lint_probe(std::int32_t count)
        : m_count{count}
    {
    }

    std::int32_t count() const
    {
        return m_count;
    }

    static std::int32_t create(std::int32_t value)
    {
        return value;
    }

private:
    std::int32_t m_count;
};

inline std::int32_t probe_twice(std::int32_t value)
{
    std::int32_t const result = style_lint_probe::create(value);
    return 2 * result;
}

}

#endif
CPP
}

replace_in_probe() {
    python3 - "$WORK/$PROBE" "$1" "$2" <<'PY'
import sys

path, old, new = sys.argv[1:]
text = open(path, encoding="utf-8").read()
if old not in text:
    sys.exit(f"the probe has no {old!r} to replace")
open(path, "w", encoding="utf-8").write(text.replace(old, new))
PY
}

plant() {
    local rule="$1" want="$2" label="$3" got
    write_probe
    if ! replace_in_probe "$4" "$5"; then
        fail "$rule: the plant could not be written"
        return
    fi
    got="$(lint_status --rule "$rule")"
    if [ "$got" -eq "$want" ]; then pass "$rule: $label (exit $got)"; else fail "$rule: $label: exit $got, wanted $want"; fi
}

before="$(git status --porcelain)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/ctrlpp-style-lint.XXXXXX")"
git ls-files -z -- lib examples scripts EXCEPTIONS.md | tar --null -T - -cf - | tar -xf - -C "$WORK"
in_copy init -q
write_probe
in_copy add -A

echo "proving the style lint can fail:"
got="$(lint_status)"
if [ "$got" -eq 0 ]; then pass "the clean probe passes every rule (exit 0)"; else fail "the clean probe: exit $got, wanted 0"; fi

plant include-order 10 "a longer standard header ahead of a shorter one" \
    $'#include <array>\n#include <cstdint>' $'#include <cstdint>\n#include <array>'
plant attributes 11 "a [[nodiscard]] on a member function" \
    '    std::int32_t count() const' '    [[nodiscard]] std::int32_t count() const'
plant return-type 12 "a trailing return type" \
    'inline std::int32_t probe_twice(std::int32_t value)' 'inline auto probe_twice(std::int32_t value) -> std::int32_t'
plant member-init 13 "a member initializer in a result struct" 'struct probe_config' 'struct probe_result'
plant member-naming 14 "a private member with a trailing underscore" 'm_count' 'count_'
plant pragma-once 15 "#pragma once" '#ifndef HPP_GUARD' $'#pragma once\n#ifndef HPP_GUARD'
plant guard 16 "a guard that does not follow the header's path" \
    'HPP_GUARD_CTRLPP_STYLE_LINT_PROBE_H' 'HPP_GUARD_STYLE_LINT_PROBE_H'
plant closing-comment 17 "a comment after the namespace's closing brace" $'}\n\n#endif' $'} // namespace ctrlpp\n\n#endif'
plant spelling 18 "British spelling in a comment" 'struct probe_config' $'// the probe\'s behaviour\nstruct probe_config'
plant discarded-result 19 "a bare call to a fallible function" \
    '    std::int32_t const result = style_lint_probe::create(value);' \
    $'    style_lint_probe::create(value);\n    std::int32_t const result = value;'
plant constant-naming 20 "a Hungarian-prefixed namespace constant" \
    $'namespace ctrlpp\n{\n' $'namespace ctrlpp\n{\n\ninline constexpr std::int32_t kProbe = 1;\n'

if [ "$(git status --porcelain)" = "$before" ]; then
    pass "the real working tree is unchanged"
else
    fail "the real working tree changed while the control ran"
fi

echo
if [ "$failures" -ne 0 ]; then
    echo "the style lint is NOT proven able to fail: $failures of the control's properties do not hold"
    exit 1
fi

echo "the style lint discriminates: the clean probe passes every rule, and each planted violation fails its"
echo "own rule with that rule's exit code"
