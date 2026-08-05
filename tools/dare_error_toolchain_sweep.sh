#!/usr/bin/env bash
#
# Driver for the discrete Riccati error-enumerator probe beside this file.
#
# It compiles and runs `dare_error_toolchain_sweep.cpp` once per configuration --
# every available compiler, by four optimization levels, by contraction off and
# on -- and prints one row per configuration over a set of bit-identical inputs.
# One input, one row per toolchain: that is what makes the table a statement
# about instruction selection rather than about the input.
#
# Each configuration is compiled and run in the same step, to its own output
# path. No binary is reused across configurations. A table whose whole content is
# the difference between configurations says nothing at all if a stale binary
# answers for one of them, and it still looks complete.
#
# ############################################################################
# # THE arm64 LEG IS ABSENT FROM THIS STATION, AND A GREEN SWEEP HERE IS NOT  #
# # EVIDENCE OF STABILITY.                                                    #
# #                                                                           #
# # Everything below is x86-64. Apple clang on arm64 exists only in continuous #
# # integration and cannot be run here. The continuous solver's twin of this   #
# # measurement produced a THIRD enumerator on arm64 that neither x86-64       #
# # compiler produced, on the same bit-identical input. A table of identical   #
# # x86-64 rows therefore establishes that these configurations agree, and     #
# # nothing whatsoever about the platform that disagreed last time.            #
# ############################################################################
#
# Usage: tools/dare_error_toolchain_sweep.sh [parallel-jobs]
#
# The compiler list is a variable at the top so an absent compiler can be dropped
# without editing the body, and the driver reports which configurations it
# actually ran rather than which it intended to.

set -u

SCRIPT_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd "${SCRIPT_DIRECTORY}/.." && pwd)"
PROBE_SOURCE="${SCRIPT_DIRECTORY}/dare_error_toolchain_sweep.cpp"

# Each entry is a compiler command followed by any flags that compiler needs to
# find a standard library it can parse. A bare clang 18 on this station selects
# gcc 16's libstdc++ and dies inside the iterator headers, which a driver that
# only checked the exit status would record as evidence about the enumerator.
COMPILER_COMMANDS=(
    "g++"
    "/usr/bin/clang++"
    "/usr/lib/llvm21/bin/clang++"
    "/usr/lib/llvm18/bin/clang++"
)
COMPILER_EXTRA_FLAGS=(
    ""
    ""
    ""
    "--gcc-install-dir=/usr/lib64/gcc/x86_64-pc-linux-gnu/15.3.0"
)

OPTIMIZATION_TAGS=("O0" "O1" "O2" "O3")

CONTRACTION_TAGS=("off" "on")
CONTRACTION_FLAGS=(
    "-ffp-contract=off"
    "-mfma -mavx2 -ffp-contract=fast"
)

# The library's own warning set for non-MSVC compilers, including both discard
# promotions, and Eigen as a system include exactly as the library's build tree
# treats it, so a diagnostic reported below comes from the project's own sources.
# The probe itself is clean under every configuration; a nonzero diagnostic count
# is the pre-existing false positive gcc raises at its first optimization level
# inside Eigen's assignment evaluator, inlined through the solver's extraction
# routine. A bare call to the solver reproduces it with this probe absent.
WARNING_FLAGS=(
    -fPIC -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
    -Wold-style-cast -Wcast-align -Woverloaded-virtual -Wnon-virtual-dtor
    -Wdouble-promotion -Wimplicit-fallthrough -Wformat=2
    -Werror=unused-result -Werror=unused-value
)
BASE_FLAGS=(-std=c++20 -fno-exceptions -fno-rtti)
INCLUDE_FLAGS=(-I "${REPOSITORY_ROOT}/lib/ctrlpp/include" -isystem /usr/include/eigen3)

PARALLEL_JOBS="${1:-3}"

WORK_DIRECTORY="${TMPDIR:-/tmp}/dare_error_toolchain_sweep"
rm -rf "${WORK_DIRECTORY}"
mkdir -p "${WORK_DIRECTORY}"

# The compiler's own reported version, reduced to a family and a release. Three
# of the four commands report themselves as clang++, and the release is what
# tells them apart.
compiler_identity()
{
    local command_path="$1"
    local first_line
    first_line="$("${command_path}" --version 2>/dev/null | head -1)"
    local release
    release="$(printf '%s' "${first_line}" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    printf '%s %s' "$(basename "${command_path}")" "${release:-unknown}"
}

anchor_field()
{
    printf '%s' "$1" | tr ' ' '\n' | grep "^$2=" | cut -d= -f2-
}

run_configuration()
{
    local index="$1"
    local command_path="$2"
    local extra_flags="$3"
    local optimization="$4"
    local contraction_flags="$5"

    local binary="${WORK_DIRECTORY}/probe_${index}"
    local diagnostics="${WORK_DIRECTORY}/diagnostics_${index}.txt"
    local result="${WORK_DIRECTORY}/result_${index}.txt"

    # shellcheck disable=SC2086
    if "${command_path}" ${extra_flags} "${BASE_FLAGS[@]}" "-${optimization}" ${contraction_flags} \
        "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" "${PROBE_SOURCE}" -o "${binary}" \
        >"${diagnostics}" 2>&1
    then
        if "${binary}" --anchor >"${result}" 2>>"${diagnostics}"; then
            return 0
        fi
        printf 'run-failed\n' >"${result}"
        return 1
    fi

    printf 'compile-failed\n' >"${result}"
    return 1
}

CONFIGURATION_COMMAND=()
CONFIGURATION_EXTRA=()
CONFIGURATION_OPTIMIZATION=()
CONFIGURATION_CONTRACTION=()
CONFIGURATION_IDENTITY=()

PLANNED_COUNT=0
SKIPPED_COUNT=0
SKIPPED_NAMES=""

for compiler_index in "${!COMPILER_COMMANDS[@]}"; do
    command_path="${COMPILER_COMMANDS[${compiler_index}]}"
    if ! command -v "${command_path}" >/dev/null 2>&1; then
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        SKIPPED_NAMES="${SKIPPED_NAMES} ${command_path}"
        continue
    fi
    identity="$(compiler_identity "${command_path}")"
    for optimization in "${OPTIMIZATION_TAGS[@]}"; do
        for contraction_index in "${!CONTRACTION_TAGS[@]}"; do
            CONFIGURATION_COMMAND+=("${command_path}")
            CONFIGURATION_EXTRA+=("${COMPILER_EXTRA_FLAGS[${compiler_index}]}")
            CONFIGURATION_OPTIMIZATION+=("${optimization}")
            CONFIGURATION_CONTRACTION+=("${contraction_index}")
            CONFIGURATION_IDENTITY+=("${identity}")
            PLANNED_COUNT=$((PLANNED_COUNT + 1))
        done
    done
done

if [ "${PLANNED_COUNT}" -eq 0 ]; then
    printf 'no compiler from the list is available; nothing was measured\n'
    exit 1
fi

# Configurations are launched in a bounded pool and every result is written to
# its own file, so the rows below are printed in a fixed order regardless of the
# order in which the builds finished.
running=0
for index in $(seq 0 $((PLANNED_COUNT - 1))); do
    run_configuration \
        "${index}" \
        "${CONFIGURATION_COMMAND[${index}]}" \
        "${CONFIGURATION_EXTRA[${index}]}" \
        "${CONFIGURATION_OPTIMIZATION[${index}]}" \
        "${CONTRACTION_FLAGS[${CONFIGURATION_CONTRACTION[${index}]}]}" &
    running=$((running + 1))
    if [ "${running}" -ge "${PARALLEL_JOBS}" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

cat <<'NOTICE'
Does the discrete Riccati solver's error enumerator depend on the input, or on
instruction selection, below its own resolution boundary?

One bit-identical set of inputs, one row per configuration. Every column is a
verdict and, where the solve was refused, the enumerator carried with it.

  accept/refuse   the two ADJACENT doubles the accept/refuse verdict changes
                  between, larger one first
  accuracy/rank   the two ADJACENT doubles the refusal changes between, from the
                  accuracy verdict to the extraction's rank verdict
  deep            fifteen binades below both, where the compared pivot ratio has
                  reached the arithmetic's own noise floor
  family-zero     the swept parameter set exactly to zero
  zero-dynamics   a state matrix and an input matrix that are both exactly zero

The first four inputs are the crossings one configuration's bisection returned,
pinned as literals rather than recomputed per configuration. A configuration
whose crossing sits elsewhere therefore reports a different verdict in those
columns instead of quietly moving the input to keep the verdict.

THE arm64 LEG IS ABSENT AND A GREEN TABLE HERE IS NOT EVIDENCE OF STABILITY.
Every row below is x86-64. Apple clang on arm64 exists only in continuous
integration and cannot be run on this station; the continuous solver's twin of
this measurement produced a third enumerator there that neither x86-64 compiler
produced. Identical rows below establish that these configurations agree, and
nothing about the platform that disagreed.

NOTICE

ROW_FORMAT='%-9s %-8s %-6s %-12s %-11s %-25s %-25s %-25s %-21s %-21s %-21s %-13s %-13s\n'
# shellcheck disable=SC2059
printf "${ROW_FORMAT}" \
    compiler release level contraction accept-side refuse-side accuracy-side rank-side \
    deep family-zero zero-dynamics accept-shift rank-shift

ROW_COUNT=0
RAN_COUNT=0
FAILED_COUNT=0
DIAGNOSTIC_LINES=0
LOWEST_BINARY=""
LOWEST_IDENTITY=""
HIGHEST_BINARY=""
HIGHEST_IDENTITY=""

ACCEPT_SIDE_VALUES=""
REFUSE_SIDE_VALUES=""
ACCURACY_SIDE_VALUES=""
RANK_SIDE_VALUES=""
DEEP_VALUES=""
FAMILY_ZERO_VALUES=""
ZERO_DYNAMICS_VALUES=""
ACCEPT_CROSSING_VALUES=""
RANK_CROSSING_VALUES=""

for index in $(seq 0 $((PLANNED_COUNT - 1))); do
    identity="${CONFIGURATION_IDENTITY[${index}]}"
    family="${identity% *}"
    release="${identity#* }"
    optimization="${CONFIGURATION_OPTIMIZATION[${index}]}"
    contraction="${CONTRACTION_TAGS[${CONFIGURATION_CONTRACTION[${index}]}]}"
    line="$(cat "${WORK_DIRECTORY}/result_${index}.txt" 2>/dev/null || printf 'missing')"
    DIAGNOSTIC_LINES=$((DIAGNOSTIC_LINES + $(wc -l <"${WORK_DIRECTORY}/diagnostics_${index}.txt" 2>/dev/null || printf 0)))

    if printf '%s' "${line}" | grep -qE '^(compile-failed|run-failed|missing)$'; then
        FAILED_COUNT=$((FAILED_COUNT + 1))
        # shellcheck disable=SC2059
        printf "${ROW_FORMAT}" \
            "${family}" "${release}" "${optimization}" "${contraction}" \
            "${line}" "${line}" "${line}" "${line}" "${line}" "${line}" "${line}" \
            "${line}" "${line}"
        ROW_COUNT=$((ROW_COUNT + 1))
        continue
    fi

    accept_side="$(anchor_field "${line}" accept_side)"
    refuse_side="$(anchor_field "${line}" refuse_side)"
    accuracy_side="$(anchor_field "${line}" accuracy_side)"
    rank_side="$(anchor_field "${line}" rank_side)"
    deep="$(anchor_field "${line}" deep)"
    family_zero="$(anchor_field "${line}" family_zero)"
    zero_dynamics="$(anchor_field "${line}" zero_dynamics)"
    accept_crossing="$(anchor_field "${line}" accept_crossing)"
    accept_shift="$(anchor_field "${line}" accept_shift)"
    rank_crossing="$(anchor_field "${line}" rank_crossing)"
    rank_shift="$(anchor_field "${line}" rank_shift)"

    # shellcheck disable=SC2059
    printf "${ROW_FORMAT}" \
        "${family}" "${release}" "${optimization}" "${contraction}" \
        "${accept_side}" "${refuse_side}" "${accuracy_side}" "${rank_side}" \
        "${deep}" "${family_zero}" "${zero_dynamics}" "${accept_shift}" "${rank_shift}"
    ROW_COUNT=$((ROW_COUNT + 1))
    RAN_COUNT=$((RAN_COUNT + 1))

    ACCEPT_SIDE_VALUES="${ACCEPT_SIDE_VALUES}${accept_side}"$'\n'
    REFUSE_SIDE_VALUES="${REFUSE_SIDE_VALUES}${refuse_side}"$'\n'
    ACCURACY_SIDE_VALUES="${ACCURACY_SIDE_VALUES}${accuracy_side}"$'\n'
    RANK_SIDE_VALUES="${RANK_SIDE_VALUES}${rank_side}"$'\n'
    DEEP_VALUES="${DEEP_VALUES}${deep}"$'\n'
    FAMILY_ZERO_VALUES="${FAMILY_ZERO_VALUES}${family_zero}"$'\n'
    ZERO_DYNAMICS_VALUES="${ZERO_DYNAMICS_VALUES}${zero_dynamics}"$'\n'
    ACCEPT_CROSSING_VALUES="${ACCEPT_CROSSING_VALUES}${accept_crossing}"$'\n'
    RANK_CROSSING_VALUES="${RANK_CROSSING_VALUES}${rank_crossing}"$'\n'

    if [ -z "${LOWEST_BINARY}" ] && [ "${optimization}" = "${OPTIMIZATION_TAGS[0]}" ] && [ "${contraction}" = "${CONTRACTION_TAGS[0]}" ]; then
        LOWEST_BINARY="${WORK_DIRECTORY}/probe_${index}"
        LOWEST_IDENTITY="${identity}, lowest optimization level, contraction ${contraction}"
    fi
    if [ -z "${HIGHEST_BINARY}" ] && [ "${optimization}" = "${OPTIMIZATION_TAGS[${#OPTIMIZATION_TAGS[@]} - 1]}" ] && [ "${contraction}" = "${CONTRACTION_TAGS[0]}" ]; then
        HIGHEST_BINARY="${WORK_DIRECTORY}/probe_${index}"
        HIGHEST_IDENTITY="${identity}, highest optimization level, contraction ${contraction}"
    fi
done

distinct_values()
{
    printf '%s' "$1" | grep -v '^$' | sort -u | tr '\n' ' '
}

distinct_count()
{
    printf '%s' "$1" | grep -v '^$' | sort -u | grep -c .
}

printf '\n'
printf 'configurations planned: %d\n' "${PLANNED_COUNT}"
printf 'configurations run: %d\n' "${RAN_COUNT}"
printf 'configurations that failed to build or run: %d\n' "${FAILED_COUNT}"
printf 'rows emitted: %d\n' "${ROW_COUNT}"
printf 'compilers skipped as unavailable: %d%s\n' "${SKIPPED_COUNT}" "${SKIPPED_NAMES}"
printf 'compiler diagnostic lines across every configuration: %d\n' "${DIAGNOSTIC_LINES}"
printf 'build logs: %s\n' "${WORK_DIRECTORY}"

printf '\ndistinct values per column, over the configurations that ran:\n'
printf '  accept-side   : %s\n' "$(distinct_values "${ACCEPT_SIDE_VALUES}")"
printf '  refuse-side   : %s\n' "$(distinct_values "${REFUSE_SIDE_VALUES}")"
printf '  accuracy-side : %s\n' "$(distinct_values "${ACCURACY_SIDE_VALUES}")"
printf '  rank-side     : %s\n' "$(distinct_values "${RANK_SIDE_VALUES}")"
printf '  deep          : %s\n' "$(distinct_values "${DEEP_VALUES}")"
printf '  family-zero   : %s\n' "$(distinct_values "${FAMILY_ZERO_VALUES}")"
printf '  zero-dynamics : %s\n' "$(distinct_values "${ZERO_DYNAMICS_VALUES}")"
printf '\ndistinct crossings each configuration located for itself:\n'
printf '  accept/refuse : %s\n' "$(distinct_values "${ACCEPT_CROSSING_VALUES}")"
printf '  accuracy/rank : %s\n' "$(distinct_values "${RANK_CROSSING_VALUES}")"

TOTAL_DISTINCT=$(( $(distinct_count "${ACCEPT_SIDE_VALUES}") \
                 + $(distinct_count "${REFUSE_SIDE_VALUES}") \
                 + $(distinct_count "${ACCURACY_SIDE_VALUES}") \
                 + $(distinct_count "${RANK_SIDE_VALUES}") \
                 + $(distinct_count "${DEEP_VALUES}") \
                 + $(distinct_count "${FAMILY_ZERO_VALUES}") \
                 + $(distinct_count "${ZERO_DYNAMICS_VALUES}") ))
COLUMN_COUNT=7

if [ "${RAN_COUNT}" -gt 0 ] && [ "${TOTAL_DISTINCT}" -eq "${COLUMN_COUNT}" ]; then
    printf '\nverdict: every column returned ONE value across every configuration that ran.\n'
    printf 'On x86-64, over these compilers and these settings, the enumerator is CONSTANT.\n'
    printf 'That is a statement about these configurations and not about arm64.\n'
else
    printf '\nverdict: at least one column returned MORE THAN ONE value across the\n'
    printf 'configurations that ran. The enumerator VARIED; the per-column lists above\n'
    printf 'name every value seen and the table names the configuration each came from.\n'
fi

if [ "${ROW_COUNT}" -ne "$((RAN_COUNT + FAILED_COUNT))" ]; then
    printf '\nrow count does not match the configurations attempted; the table is incomplete\n'
    exit 1
fi

if [ "${FAILED_COUNT}" -ne 0 ]; then
    printf '\nat least one configuration did not build or run; see the build logs above\n'
    exit 1
fi

if [ -n "${HIGHEST_BINARY}" ] && [ -n "${LOWEST_BINARY}" ]; then
    "${HIGHEST_BINARY}" >"${WORK_DIRECTORY}/full_highest.txt"
    "${LOWEST_BINARY}" >"${WORK_DIRECTORY}/full_lowest.txt"
    if cmp -s "${WORK_DIRECTORY}/full_highest.txt" "${WORK_DIRECTORY}/full_lowest.txt"; then
        identical="yes"
    else
        identical="no"
    fi

    printf '\n== full sweep, from one configuration of the table above ==\n'
    printf 'reference configuration: %s\n' "${HIGHEST_IDENTITY}"
    printf 'byte-identical to the same sweep from %s: %s\n\n' "${LOWEST_IDENTITY}" "${identical}"
    cat "${WORK_DIRECTORY}/full_highest.txt"
fi
