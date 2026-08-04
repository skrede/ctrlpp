#!/usr/bin/env bash
#
# Source purity gate for the shipped ctrlpp surface.
#
# Seven rules, checked from one place so a single run tells a developer
# everything that must be fixed. Every failing rule is reported before the
# script exits; it never stops at the first.
#
#   Rule 1  No unwrap accessor. Scans lib/ and examples/. Bans the accessor for
#           the optional type and for the library's own fallible-result backport
#           together: separating the two in a text search is the fiddly part,
#           and an explicit branch is the correct style for both. Examples are
#           documentation users copy, which is why they are in scope. The
#           backport header is exempt at FILE level, because it defines the
#           accessor and documents it.
#
#   Rule 2  No throw statement. Scans lib/ and examples/. The backport header is
#           exempt at FILE level, because it holds two throw statements -- the
#           gated one its contract names and a bare rethrow beside it. A
#           line-level exemption would miss the second.
#
#   Rule 3  No run-time type identification. Scans lib/ and examples/. Already
#           at zero; the rule makes that permanent.
#
#   Rule 4a The class-level discard-warning attribute is PRESENT on every
#           declaration in the manifest below. Scans lib/. This is a count
#           assertion, not a text search: an attribute deleted from a header
#           makes it fail, which a search for the attribute could never detect.
#
#   Rule 4b No SITE-level discard-warning attribute exists outside the
#           manifest's declarations. Scans lib/, tests/, examples/, benchmarks/
#           and validation/. This is the scope fence. The attribute is earned
#           only where discarding the result silently breaks correctness -- a
#           dropped failure, or a dropped substitution report, inside a control
#           loop. It is forbidden on trivial accessors, getters, size and empty
#           queries, and state-health predicates. Widening it beyond the failure
#           channel and the disposition channel is a decision for the user, not
#           a decision this gate may make on its own.
#
#   Rule 4c No bare-statement call to a manifest-listed fallible function.
#           Scans lib/ and examples/. THIS RULE IS A RATCHET, NOT THE
#           ENFORCEMENT. The enforcement is the compiler: the discard warning is
#           promoted to an error on ctrlpp's own targets, so a discarded
#           fallible return does not build. The configuration-time compile
#           canary in tests/CMakeLists.txt proves this with the active compiler.
#           Rule 4c reaches only what the compiler does not see -- a call site
#           in a target that neither build tree compiles.
#
#   Rule 5  The exception-mode macro appears in no file under lib/ but two: the
#           fallible-result backport header, and the configuration header that
#           defines the macro.
#
#   Rule 6  The flags that promote the discard warning to an error remain
#           declared in the top-level build file for both compiler families.
#           This is a cheap cross-platform diagnostic, not proof of compiler
#           enforcement. The positive-and-negative compile canary configured
#           from tests/CMakeLists.txt supplies that proof for the active
#           compiler.
#
#   Rule 7  The static-analysis check that re-adds site-level discard-warning
#           attributes on trivial returns stays disabled. Enabling it would
#           reintroduce, mechanically, exactly what rule 4b forbids.
#
#   Rule 8  No planning-artifact identifier appears in shipped source,
#           documentation, examples, benchmarks or tests. Decision-record keys,
#           seed keys, requirement and task keys, planning-document filenames
#           and paths into the planning directory are all forbidden by the
#           project's own rules: a reader of this repository must not need
#           access to a planning system to understand a comment, and the
#           artifact a key points at is not shipped with the code. Scans every
#           file under the roots, not only the C++ ones, because one such key
#           was found in a build file.
#
# Usage:
#   scripts/source_purity_check.sh [scan-root]
#   scripts/source_purity_check.sh --self-test
#   PURITY_SCAN_ROOT=some/root scripts/source_purity_check.sh
#
# Environment:
#   PURITY_SCAN_ROOT   directory the rules scan (default: the repository root)
#
# --- What this gate does NOT establish ----------------------------------------
#
# Stated here rather than left implicit, because a gate that hides its blind
# spots is the same defect in a different place.
#
#   1. These rules are text searches. They are a ratchet against regression, not
#      a proof that no reachable path throws. Nothing here is evidence about
#      reachability, and a passing run must not be presented as one.
#
#   2. Rule 2 drops lines whose first non-space character opens a comment. A
#      real throw placed after code on a line that ALSO opens a trailing comment
#      would be dropped with it. That shape does not occur today and is not
#      worth a C++ parser, but it is a known hole. The tempting alternative --
#      anchoring on the first character not being a slash -- is deliberately NOT
#      used: it silently drops any line containing a division operator, and a
#      gate that passes for the wrong reason is worse than no gate.
#
#   3. This script necessarily contains the patterns it searches for, and sits
#      under scripts/, which is outside every scan root. That exclusion is
#      deliberate, not accidental: a gate that scanned its own implementation
#      would flag itself and could never pass.
#
#   4. Rule 4a checks the manifest's entries and nothing else. A new failure
#      channel or disposition type added without a manifest row is not checked.
#      The manifest is maintained by hand, and that is the cost of a list-driven
#      rule. Two report types in the model-predictive and moving-horizon modules
#      are disposition-shaped and carry no class-level attribute today; they are
#      deliberately absent from the manifest, which records what shipped rather
#      than what a reader might expect, and annotating them is a user decision.
#
#   5. Rule 4c is a ratchet behind the compiler, not the enforcement. The
#      configuration-time compile canary proves that a valid translation unit
#      succeeds without promotion and fails when the active compiler's
#      discarded-result warning is promoted. Rule 6 only confirms that both
#      compiler-family spellings remain declared; a reader must not take that
#      text search for the compiler's work.
#
#   6. Rule 4c is also name-based, so it can only carry names that are fallible
#      on EVERY declaration in the library. Names shared with an infallible
#      surface are excluded and listed beside the manifest with the reason. For
#      those, the compiler is the whole enforcement. Rule 4c further assumes a
#      statement fits on one line: a call whose argument list wraps is not seen,
#      and a call whose result is bound on the line above is suppressed by a
#      continuation test rather than parsed.
#
#   7. The roots are lib/, tests/, examples/, benchmarks/ and validation/, plus
#      docs/ for rule 8 ONLY. The documentation tree is deliberately not a root
#      for the attribute rules, and that choice is recorded rather than
#      inherited: under the current attribute policy a documentation page that
#      shows the class-level attribute on a quoted signature is showing correct
#      documentation, so scanning the documentation tree for the attribute would
#      be actively wrong. Rule 8 is the opposite case -- the project's rule names
#      documentation explicitly -- so it scans there. Rules 1 to 3 do not reach
#      the documentation tree; a documentation page that unwraps teaches the
#      idiom rule 1 bans, and closing that is a separate, wider pass with this
#      script's own path excluded by NAME rather than by extension -- extension
#      is precisely the exclusion that lets documentation pages through.
#
#   8. Rule 8 matches identifier FORMS, never the English words around them.
#      That is a measurement, not a preference: a pattern keyed on the word
#      "phase" reported 41 hits of which 20 were motion-profile phases, a stage
#      of a reference linear-algebra routine, and tutorial headings. It carries
#      two consequences. A tagged key needs TWO digits to match, because one
#      digit is indistinguishable from template-dimension arithmetic such as
#      NX-1; single-digit keys are matched only for the specific families known
#      to use them. And two designations that match the shape are excluded by
#      NAME, listed beside the rule: they are an external floating-point
#      standard and a board form factor, neither of which is a planning key.
#      That exclusion is exercised by the real tree rather than by a fixture,
#      since both designations are present in it today.
#
# This is a plain, repeatable local script. It authors no continuous-integration
# configuration and no CMake toolchain file.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

# --- Arguments ----------------------------------------------------------------

self_test=0
scan_root_arg=""

for arg in "$@"; do
    case "${arg}" in
        --self-test)
            self_test=1
            ;;
        -h|--help)
            sed -n '2,70p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        -*)
            echo "ERROR: unknown option '${arg}'." >&2
            echo "  Usage: scripts/source_purity_check.sh [scan-root]" >&2
            echo "         scripts/source_purity_check.sh --self-test" >&2
            exit 2
            ;;
        *)
            scan_root_arg="${arg}"
            ;;
    esac
done

fixture_root="${repo_root}/scripts/purity_fixtures"
self_test_root=""

if [ "${self_test}" -eq 1 ]; then
    if [ -n "${scan_root_arg}" ]; then
        echo "ERROR: --self-test scans its own fixture directory and takes no root." >&2
        exit 2
    fi

    self_test_root="$(mktemp -d "${TMPDIR:-/tmp}/ctrlpp-purity-self-test.XXXXXX")"
    cleanup_self_test()
    {
        rm -rf -- "${self_test_root}"
    }
    trap cleanup_self_test EXIT

    cp -R "${fixture_root}/." "${self_test_root}/"

    decision_record_key="D""-07"
    versioned_branch_key="mile""stone/v7.8.9"
    printf '// Deliberate self-test sample: %s\n' "${decision_record_key}" \
        > "${self_test_root}/lib/generated_decision_record_identifier.h"
    printf '// Deliberate self-test sample: %s\n' "${versioned_branch_key}" \
        > "${self_test_root}/lib/generated_versioned_branch_identifier.h"

    scan_root="${self_test_root}"
else
    scan_root="${scan_root_arg:-${PURITY_SCAN_ROOT:-${repo_root}}}"
fi

if [ ! -d "${scan_root}" ]; then
    echo "ERROR: scan root '${scan_root}' does not exist." >&2
    echo "  Pass a directory holding lib/, tests/, examples/, benchmarks/ or" >&2
    echo "  validation/, or run with no argument to scan the repository root." >&2
    exit 2
fi

scan_root="$(cd "${scan_root}" && pwd)"
cd "${scan_root}"

# --- The manifest -------------------------------------------------------------
#
# One row per declaration that MUST carry the class-level discard-warning
# attribute. Each row is a path relative to the scan root, a vertical bar, and
# the exact declaration text the attribute sits on with surrounding whitespace
# trimmed.
#
# Rule 4a asserts each row is present exactly once. Rule 4b treats these rows,
# and only these rows, as the permitted places for the attribute to appear.
#
# Built by hand from what each contributing change shipped, and verified against
# the headers rather than transcribed from a forecast. It is deliberately NOT
# built by searching for the attribute: a search finds what IS annotated, which
# cannot detect an annotation that was deleted, and detecting that deletion is
# the entire point of rule 4a.

manifest=(
    # The fallible result type callers branch on. One annotation here covers
    # every fallible return in the library, present and future.
    "lib/ctrlpp/include/ctrlpp/detail/expected.h|class [[nodiscard]] expected"
    # Its void specialization -- a partial specialization does not inherit the
    # primary template's attribute, which makes this the row most likely to be
    # forgotten and the reason 4a is a count rather than a presence search.
    "lib/ctrlpp/include/ctrlpp/detail/expected.h|class [[nodiscard]] expected<void, E>"
    # What the online planners realized against what they were commanded, when
    # the planner brakes and replans instead. Discarding it in a control loop
    # loses the substitution silently.
    "lib/ctrlpp/include/ctrlpp/trajectory/online_planner_diagnostics.h|struct [[nodiscard]] online_planner_diagnostics"
    # What the trapezoidal profile realized against what it was commanded, when
    # boundary velocities force the acceleration up. Same reason.
    "lib/ctrlpp/include/ctrlpp/trajectory/trapezoidal_trajectory.h|struct [[nodiscard]] trapezoidal_disposition"
)

# The self-test's twin manifest: same shape, paths under the fixture root. Rule
# 4a is a presence assertion, so it needs a manifest that resolves against the
# fixture in order to be exercised at all. Excluding it from the self-test would
# leave the one rule most likely to be skipped -- on the argument that a
# presence check is awkward to violate -- with no evidence that it can fail.
#
# Two rows, both in one file, and the pairing is load-bearing. The first is
# SATISFIED and the second is not, so the self-test's expected count for rule 4a
# is one. Deleting the fixture file then breaks BOTH rows and the count goes to
# two, which is what lets the self-test tell a deleted fixture apart from a
# deleted annotation. With only the unsatisfied row, a missing file and an
# absent attribute produce the identical signal and the fixture could be deleted
# with the self-test still green. The satisfied row also exercises rule 4a's
# positive branch and rule 4b's manifest exemption, neither of which any other
# fixture reaches.
fixture_manifest=(
    "lib/fixture_missing_class_annotation.h|struct [[nodiscard]] fixture_present_disposition"
    "lib/fixture_missing_class_annotation.h|struct [[nodiscard]] fixture_disposition"
)

if [ "${self_test}" -eq 1 ]; then
    manifest=("${fixture_manifest[@]}")
fi

# --- Rule 4c's callable list --------------------------------------------------
#
# Function and member names whose return is fallible on EVERY declaration in the
# library. A bare-statement call to one of these drops a typed failure.
#
# Deliberately EXCLUDED, because each name is shared with an infallible surface
# and a name-based search cannot tell them apart from a call site:
#
#   update    fallible on the six estimators and the complementary filter;
#             returns nothing on the particle filter, on both online planners,
#             on the null observer, and inside the observer concept's
#             requires-expression
#   compute   fallible on the proportional-integral-derivative controller;
#             returns a control vector on the linear-quadratic gain
#   evaluate  fallible on the two adaptive controllers; returns a trajectory
#             point across the whole trajectory family
#   solve     fallible on the predictive controllers; returns a solver result
#             aggregate on the solver adapters
#   normalize fallible on the rotation-group helper; the linear-algebra
#             dependency's in-place member of the same name returns nothing
#
# For those five the compiler is the whole enforcement. Adding them here would
# require an exemption for every correct call to an infallible surface, and a
# rule whose exemption list grows as correct code is written is not a ratchet.

fallible_call_names="create|setup|place_observer|place|lqr_gain_continuous|lqr_gain|partition_lqi_gain|lqi_gain|validate_biquad_design|can_rescale_to|rescale_to|dare|care"

# --- Rule 8's identifier forms ------------------------------------------------
#
# Keys produced by a planning system, matched by SHAPE. In order: decision
# records, audit findings, success criteria, seed documents, threat entries,
# the general requirement-or-task key (an uppercase tag, a dash, and at least
# TWO digits), planning-document filenames, and any path into the planning
# directory.
#
# The two-digit floor on the general form is what keeps template-dimension
# arithmetic out: NX-1, NU-1 and NB-1 all have the shape but one digit. The
# families that legitimately use a single digit are spelled out individually
# instead.
planning_identifier_pattern='(^|[^A-Za-z0-9_-])(D-[0-9]+|DL-[0-9]+|SC-?[0-9]+|SEED-[0-9]+|T-[0-9]+-[0-9]+-[0-9]+|[A-Z][A-Z0-9]{1,9}-[0-9]{2,}|milestone/[A-Za-z0-9][A-Za-z0-9._-]*|(RESEARCH|PLAN|SUMMARY|CONTEXT|ROADMAP|PATTERNS)\.md|\.planning/)'

# Designations that match the general shape and are not planning keys, excluded
# by name rather than by widening or narrowing the pattern:
#   IEEE-754    the binary floating-point standard
#   NUCLEO-144  the board form factor of the embedded example's target
#   SC2086      a shellcheck directive code, which has the success-criteria
#               shape and none of the meaning
#   .planning/  as a YAML list entry under a workflow's paths-ignore, which is
#               the one legitimate mention of that directory in shipped
#               material: it tells the runner to skip planning-only pushes
#               rather than pointing a reader at an artifact they cannot read.
#               Anchored to the list-entry form, so a .planning/ path anywhere
#               else in a workflow is still a violation.
# All four are present in the tree today, so a real run is what proves these
# exclusions still work; no fixture is needed for them and none would be honest.
planning_identifier_exclusions="(IEEE|NUCLEO)-[0-9]|shellcheck[[:space:]]+disable=SC[0-9]|:[0-9]+:[[:space:]]*-[[:space:]]*['\"]?\\.planning/"

# --- Shared helpers -----------------------------------------------------------

# Drops lines whose first non-space character opens a comment: a double slash, a
# documentation triple slash, a block-continuation star, or a block opener.
comment_filter()
{
    grep -vE ':[0-9]+:[[:space:]]*(//|\*|/\*)' || true
}

scan_dirs=()

# Collects those of the named directories that exist under the scan root, so a
# fixture root holding only one of them does not make a rule error out.
select_scan_dirs()
{
    local d
    scan_dirs=()
    for d in "$@"; do
        if [ -d "${d}" ]; then
            scan_dirs+=("${d}")
        fi
    done
    return 0
}

# BUILD TREES ARE NOT SCANNED, and the reason is what this gate is for.
#
# Every root the rules select is a directory of TRACKED material. A build tree
# sitting inside one of them is generated output, and when a dependency is
# fetched rather than found installed -- which is what the benchmark project's
# own CTRLPP_FETCH_BENCHMARK_DEPS option does, and the only way to build the
# benches on a machine without the dependency installed -- that output contains
# a third-party library's sources and its rendered documentation. Reading it can
# only produce findings about someone else's code. This gate exists to catch a
# false claim in material this repository ships; a report about a vendored
# library's analytics tag or its own nodiscard macro is the false red the gate is
# supposed to prevent, not one it should manufacture.
#
# The patterns mirror the .gitignore entries that create these directories, so
# what the gate declines to read is exactly what the repository declines to
# track. Verified: no tracked path under any scanned root matches them.
#
# Applied to EVERY recursive search rather than to the two rules that scan
# `benchmarks` today, so that a rule which later adds a root carrying a build
# tree inherits the exclusion instead of rediscovering this failure.
scan_exclusions=(--exclude-dir='build*' --exclude-dir='cmake-build-*')

# Counts the lines of a file whose whitespace-trimmed text equals a given
# string. An exact match, so a mention of the same text inside a comment does
# not satisfy the assertion.
count_exact_lines()
{
    local file="$1" decl="$2"
    if [ ! -f "${file}" ]; then
        printf '0\n'
        return 0
    fi
    awk -v d="${decl}" '
        {
            l = $0
            sub(/^[[:blank:]]+/, "", l)
            sub(/[[:blank:]]+$/, "", l)
            if (l == d) { n++ }
        }
        END { print n + 0 }
    ' "${file}"
}

trim()
{
    printf '%s' "$1" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

# --- The rules ----------------------------------------------------------------
#
# Each rule writes one line per violation to standard output and nothing else.
# The driver counts the lines.

rule_1_unwrap()
{
    select_scan_dirs lib examples
    [ ${#scan_dirs[@]} -eq 0 ] && return 0
    grep -rn --include='*.h' --include='*.hpp' --include='*.cpp' -E '\.value\(\)' \
        "${scan_exclusions[@]}" "${scan_dirs[@]}" 2>/dev/null \
        | comment_filter \
        | grep -v '/detail/expected\.h:' \
        || true
    return 0
}

rule_2_throw()
{
    select_scan_dirs lib examples
    [ ${#scan_dirs[@]} -eq 0 ] && return 0
    grep -rn --include='*.h' --include='*.hpp' --include='*.cpp' \
        -E '(^|[^[:alnum:]_])throw([[:space:]]|;|$)' \
        "${scan_exclusions[@]}" "${scan_dirs[@]}" 2>/dev/null \
        | comment_filter \
        | grep -v '/detail/expected\.h:' \
        || true
    return 0
}

rule_3_rtti()
{
    select_scan_dirs lib examples
    [ ${#scan_dirs[@]} -eq 0 ] && return 0
    grep -rn --include='*.h' --include='*.hpp' --include='*.cpp' \
        -E '(^|[^[:alnum:]_])(dynamic_cast|typeid)[[:space:]]*[<(]' \
        "${scan_exclusions[@]}" "${scan_dirs[@]}" 2>/dev/null \
        | comment_filter \
        || true
    return 0
}

rule_4a_class_attribute_present()
{
    local entry mpath mdecl found
    for entry in ${manifest[@]+"${manifest[@]}"}; do
        mpath="${entry%%|*}"
        mdecl="${entry#*|}"
        if [ ! -f "${mpath}" ]; then
            printf '%s: manifest file is missing, so the declaration "%s" cannot be verified\n' \
                "${mpath}" "${mdecl}"
            continue
        fi
        found="$(count_exact_lines "${mpath}" "${mdecl}")"
        if [ "${found}" != "1" ]; then
            printf '%s: expected exactly 1 declaration reading "%s", found %s\n' \
                "${mpath}" "${mdecl}" "${found}"
        fi
    done
    return 0
}

rule_4b_no_site_attribute()
{
    local hits hit file rest lineno content trimmed entry mpath mdecl exempt
    select_scan_dirs lib tests examples benchmarks validation
    [ ${#scan_dirs[@]} -eq 0 ] && return 0
    hits="$(grep -rn --include='*.h' --include='*.hpp' --include='*.cpp' \
        -F '[[nodiscard]]' "${scan_exclusions[@]}" "${scan_dirs[@]}" 2>/dev/null | comment_filter || true)"
    while IFS= read -r hit; do
        [ -n "${hit}" ] || continue
        file="${hit%%:*}"
        rest="${hit#*:}"
        lineno="${rest%%:*}"
        content="${rest#*:}"
        trimmed="$(trim "${content}")"
        exempt=0
        for entry in ${manifest[@]+"${manifest[@]}"}; do
            mpath="${entry%%|*}"
            mdecl="${entry#*|}"
            if [ "${file}" = "${mpath}" ] && [ "${trimmed}" = "${mdecl}" ]; then
                exempt=1
                break
            fi
        done
        if [ "${exempt}" -eq 0 ]; then
            printf '%s:%s: site-level attribute outside the manifest: %s\n' \
                "${file}" "${lineno}" "${trimmed}"
        fi
    done <<EOF
${hits}
EOF
    return 0
}

rule_4c_discarded_fallible()
{
    select_scan_dirs lib examples
    [ ${#scan_dirs[@]} -eq 0 ] && return 0
    find "${scan_dirs[@]}" -type f \
        \( -name '*.h' -o -name '*.hpp' -o -name '*.cpp' \) \
        -exec awk -v names="${fallible_call_names}" '
        FNR == 1 { prev = "" }
        {
            line = $0
            s = line
            sub(/^[[:blank:]]+/, "", s)
            sub(/[[:blank:]]+$/, "", s)
            if (s == "") { next }
            # A comment line is neither a statement nor a predecessor of one.
            if (s ~ /^(\/\/|\*|\/\*)/) { next }
            call_re = "^([A-Za-z_][A-Za-z0-9_]*(\\.|->|::))*(" names ")[[:blank:]]*\\("
            is_call = (s ~ call_re)
            # A complete statement ends the line, so the argument list closed.
            ends_statement = (s ~ /;[[:blank:]]*(\/\/.*)?$/)
            # A trailing return type marks a declaration, never a call.
            is_declaration = (s ~ /\)[[:blank:]]*(const[[:blank:]]*)?(noexcept[[:blank:]]*)?->/)
            # A line whose predecessor left an expression open is a
            # continuation: the result is bound above, so nothing is discarded.
            is_continuation = (prev ~ /[=,(]$/ || prev ~ /(^|[^A-Za-z0-9_])return$/)
            if (is_call && ends_statement && !is_declaration && !is_continuation) {
                printf "%s:%d: fallible result discarded: %s\n", FILENAME, FNR, s
            }
            prev = s
        }
    ' {} +
    return 0
}

rule_5_exception_macro()
{
    local file
    [ -d lib ] || return 0
    for file in $(grep -rl 'CTRLPP_HAS_EXCEPTIONS' lib 2>/dev/null || true); do
        case "${file}" in
            */detail/expected.h|*/config.h)
                ;;
            *)
                printf '%s: exception-mode macro outside the backport header and the configuration header that defines it\n' \
                    "${file}"
                ;;
        esac
    done
    return 0
}

rule_6_promotion_flag()
{
    local build_file="CMakeLists.txt" token
    if [ ! -f "${build_file}" ]; then
        printf '%s: missing, so the discard-warning promotion cannot be verified\n' "${build_file}"
        return 0
    fi
    for token in '-Werror=unused-result' '/we4834'; do
        if ! grep -qF -- "${token}" "${build_file}"; then
            printf '%s: the discard-warning promotion "%s" is absent\n' \
                "${build_file}" "${token}"
        fi
    done
    return 0
}

rule_7_analysis_check_disabled()
{
    local config_file=".clang-tidy"
    if [ ! -f "${config_file}" ]; then
        printf '%s: missing, so the attribute-inserting check cannot be verified as disabled\n' \
            "${config_file}"
        return 0
    fi
    grep -nE '^[[:space:]]*modernize-use-nodiscard' "${config_file}" \
        | sed "s|^|${config_file}:|; s|\$| <- the attribute-inserting check is enabled|" \
        || true
    return 0
}

rule_8_planning_identifiers()
{
    # .github is scanned for the same reason the roots below are, and it was
    # added because two seed keys were living in the fuzz workflow's waiver list
    # while this rule passed. Continuous-integration configuration is material a
    # reader of this repository reads, and a waiver that names a key rather than
    # the defect it waives tells that reader nothing they can act on.
    select_scan_dirs lib tests examples benchmarks validation docs .github
    [ ${#scan_dirs[@]} -eq 0 ] && return 0
    # Every file, not only the C++ ones: one such key was found in a build file,
    # which an extension filter would have missed. Binary files are skipped,
    # because a key inside one carries nothing to a reader.
    #
    # There is no comment filter here, deliberately. Every occurrence this rule
    # was written against was inside a comment; a comment is the place these
    # keys live, not an exemption from them.
    grep -rInE "${planning_identifier_pattern}" \
        "${scan_exclusions[@]}" "${scan_dirs[@]}" 2>/dev/null \
        | grep -vE "${planning_identifier_exclusions}" \
        || true
    return 0
}

# --- Driver -------------------------------------------------------------------

rule_ids=()
rule_counts=()
rule_expected_counts=()
total_violations=0

run_rule()
{
    local id="$1" title="$2" fn="$3"
    local out count

    echo "=== Rule ${id}: ${title} ==="
    out="$("${fn}")"
    count=0
    if [ -n "${out}" ]; then
        count="$(printf '%s\n' "${out}" | wc -l | tr -d '[:space:]')"
        printf '%s\n' "${out}" | sed 's/^/  /'
    fi

    local expected_count=1
    if [ "${id}" = "8" ]; then
        expected_count=2
    fi

    rule_ids+=("${id}")
    rule_counts+=("${count}")
    rule_expected_counts+=("${expected_count}")
    total_violations=$((total_violations + count))

    if [ "${self_test}" -eq 1 ]; then
        echo "Rule ${id}: ${count} violation(s) detected; expected ${expected_count}"
    elif [ "${count}" -eq 0 ]; then
        echo "Rule ${id} PASS"
    else
        echo "Rule ${id} FAIL (${count} violation(s))"
    fi
    echo
}

echo "Scan root: ${scan_root}"
if [ "${self_test}" -eq 1 ]; then
    echo "Mode:      self-test against the deliberately-violating fixture"
fi
echo

run_rule "1"  "no unwrap accessor (lib, examples)"                       rule_1_unwrap
run_rule "2"  "no throw statement (lib, examples)"                       rule_2_throw
run_rule "3"  "no run-time type identification (lib, examples)"          rule_3_rtti
run_rule "4a" "class-level discard attribute present (lib)"              rule_4a_class_attribute_present
run_rule "4b" "no site-level discard attribute outside the manifest"     rule_4b_no_site_attribute
run_rule "4c" "no bare-statement call to a fallible function"            rule_4c_discarded_fallible
run_rule "5"  "exception-mode macro confined to two files (lib)"         rule_5_exception_macro
run_rule "6"  "discard-error flags declared for both compiler families"  rule_6_promotion_flag
run_rule "7"  "attribute-inserting analysis check stays disabled"        rule_7_analysis_check_disabled
run_rule "8"  "no planning-artifact identifier in shipped material"     rule_8_planning_identifiers

if [ "${self_test}" -eq 1 ]; then
    # The exact-count assertion is what makes this a self-test rather than a
    # check that something went wrong: a script that crashed before scanning
    # would also exit nonzero, and only the count distinguishes the two.
    self_test_failed=0
    expected_total=0
    i=0
    while [ "${i}" -lt "${#rule_ids[@]}" ]; do
        expected_total=$((expected_total + rule_expected_counts[i]))
        if [ "${rule_counts[${i}]}" -ne "${rule_expected_counts[${i}]}" ]; then
            echo "SELF-TEST FAIL: rule ${rule_ids[${i}]} reported ${rule_counts[${i}]} violation(s), expected ${rule_expected_counts[${i}]}." >&2
            if [ "${rule_counts[${i}]}" -eq 0 ]; then
                echo "  Its fixture produced nothing: either the fixture file is gone or the rule stopped detecting." >&2
            else
                echo "  Its fixture produced more than one hit: either a filter the rule depends on was removed," >&2
                echo "  or a fixture file it reads is gone. The reported lines above say which." >&2
            fi
            self_test_failed=1
        fi
        i=$((i + 1))
    done
    if [ "${self_test_failed}" -ne 0 ]; then
        exit 1
    fi
    echo "Self-test PASS: ${total_violations} violation(s) detected, matching all ${#rule_ids[@]} rule expectations."
    exit 0
fi

if [ "${total_violations}" -gt 0 ]; then
    echo "${total_violations} violation(s) across $(printf '%s\n' ${rule_counts[@]+"${rule_counts[@]}"} | grep -vc '^0$' || true) rule(s)."
    echo "  Every rule above reports independently; fix each FAIL section."
    exit 1
fi

echo "All rules PASS across ${#rule_ids[@]} rules."
