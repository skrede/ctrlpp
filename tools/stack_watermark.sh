#!/usr/bin/env bash
#
# Driver for the stack watermark instrument beside this file.
#
# It compiles and runs `stack_watermark.cpp` once per state dimension, under the
# provenance the published Riccati stack figures were taken at, and prints two
# measurements side by side with the published figures next to them.
#
# ##########################################################################
# # THE TWO MEASUREMENTS ARE NOT INTERCHANGEABLE AND NEITHER SUBSTITUTES   #
# # FOR THE OTHER.                                                         #
# #                                                                        #
# # The compiler's per-function frame report attributes NOTHING TO CALLEES. #
# # A sum over named functions is therefore a LOWER BOUND on what a task    #
# # stack must hold: it excludes every Eigen and libm frame beneath the     #
# # named ones. The runtime watermark covers the whole chain and is the     #
# # number a task stack must actually cover. Reporting the frame figure as  #
# # a row's stack cost understates it, and a heap-only or frame-only        #
# # statement is blind to what the stack has to hold.                       #
# ##########################################################################
#
# Each dimension is compiled and run in the same step, to its own output path.
# No binary is reused across dimensions. A table whose whole content is the
# difference between dimensions says nothing at all if a stale binary answers for
# one of them, and it still looks complete.
#
# The comparison against the published figures is EXACT and carries no
# tolerance. A harness adjusted until the figures agree is not a control, so a
# mismatch is printed as a mismatch and the driver exits nonzero. So does a
# nonzero harness floor, which means the figure beside it is not attributable to
# the call under measurement.
#
# Usage: tools/stack_watermark.sh [--frames] [--input-dimension N] [--jobs N]
#
#   --frames             also take the per-function frame measurement, which
#                        costs a second compile per dimension
#   --input-dimension N  the corpus input dimension (default 1)
#   --jobs N             concurrent compiles (default 2)

set -u

SCRIPT_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd "${SCRIPT_DIRECTORY}/.." && pwd)"
INSTRUMENT_SOURCE="${SCRIPT_DIRECTORY}/stack_watermark.cpp"

# The published provenance, reproduced flag for flag: g++, C++20, -O2, and NO
# architecture flag so the driver default applies. Changing any of these makes
# the rows below incomparable with the figures they are checked against, which is
# the reason they are written here rather than passed in.
COMPILER_COMMAND="g++"
BASE_FLAGS=(-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread)

# The library's own warning set for non-MSVC compilers, including both discard
# promotions, and Eigen as a system include exactly as the library's build tree
# treats it, so a diagnostic reported below comes from the instrument's own
# source.
WARNING_FLAGS=(
    -fPIC -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
    -Wold-style-cast -Wcast-align -Woverloaded-virtual -Wnon-virtual-dtor
    -Wdouble-promotion -Wimplicit-fallthrough -Wformat=2
    -Werror=unused-result -Werror=unused-value
)
INCLUDE_FLAGS=(-I "${REPOSITORY_ROOT}/lib/ctrlpp/include" -isystem /usr/include/eigen3)

STATE_DIMENSIONS=(2 4 6 8)

# The whole-chain runtime peaks and the acceptance-check frames the real-time
# safety matrix publishes for the discrete Riccati row, keyed by state dimension.
declare -A PUBLISHED_PEAK=([2]=5352 [4]=11688 [6]=23144 [8]=42760)
declare -A PUBLISHED_FRAME=([2]=864 [4]=2400 [6]=6432 [8]=15056)

# The published frame column carries THIS function's frame, and it is not in
# general the deepest frame in the chain: at two states a Schur reordering
# routine has a larger one. The driver prints both rather than letting either
# stand for the other.
PUBLISHED_FRAME_FUNCTION="estimate_riccati_forward_error"

INPUT_DIMENSION=1
PARALLEL_JOBS=2
WITH_FRAMES=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --frames)
            WITH_FRAMES=1
            shift
            ;;
        --input-dimension)
            INPUT_DIMENSION="$2"
            shift 2
            ;;
        --jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        *)
            printf 'unrecognized argument: %s\n' "$1" >&2
            exit 2
            ;;
    esac
done

if ! command -v "${COMPILER_COMMAND}" >/dev/null 2>&1; then
    printf '%s is not available; nothing was measured\n' "${COMPILER_COMMAND}" >&2
    exit 1
fi

# Named for the sweep rather than for the instrument, so that the directory does
# not collide with a single-configuration binary a developer builds by hand at
# the obvious path.
WORK_DIRECTORY="${TMPDIR:-/tmp}/stack_watermark_sweep"
rm -rf "${WORK_DIRECTORY}"
mkdir -p "${WORK_DIRECTORY}"

COMPILER_VERSION="$("${COMPILER_COMMAND}" --version 2>/dev/null | head -1)"
HOST_ARCHITECTURE="$(uname -m)"
HOST_SYSTEM="$(uname -s)"

# The thread stack the instrument requests, read out of the instrument itself so
# the provenance block cannot drift away from what is actually measured.
THREAD_STACK_MIB="$(grep -oE 'thread_stack_bytes = std::size_t\{[0-9]+\}' "${INSTRUMENT_SOURCE}" |
    grep -oE '[0-9]+' | tail -1)"

measure_dimension()
{
    local state_dimension="$1"
    local binary="${WORK_DIRECTORY}/instrument_${state_dimension}"
    local diagnostics="${WORK_DIRECTORY}/diagnostics_${state_dimension}.txt"
    local result="${WORK_DIRECTORY}/watermark_${state_dimension}.txt"

    if "${COMPILER_COMMAND}" "${BASE_FLAGS[@]}" "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" \
        "-DWATERMARK_STATE_DIMENSION=${state_dimension}" \
        "-DWATERMARK_INPUT_DIMENSION=${INPUT_DIMENSION}" \
        "${INSTRUMENT_SOURCE}" -o "${binary}" >"${diagnostics}" 2>&1
    then
        if ! "${binary}" >"${result}" 2>>"${diagnostics}"; then
            printf 'run-failed\n' >"${result}"
        fi
    else
        printf 'compile-failed\n' >"${result}"
    fi
}

# The per-function frame report. The compiler writes one line per EMITTED
# function, so a function the optimizer inlined has no line of its own and its
# frame is merged into its caller's. That is why the instrument's own
# translation unit is excluded from the "deepest" search below: the solve is
# inlined into the instrument's call wrapper at the larger dimensions, and
# reporting the wrapper would attribute the library's frame to the harness.
measure_frames()
{
    local state_dimension="$1"
    local object="${WORK_DIRECTORY}/frames_${state_dimension}.o"
    local report="${WORK_DIRECTORY}/frames_${state_dimension}.su"
    local diagnostics="${WORK_DIRECTORY}/frame_diagnostics_${state_dimension}.txt"
    local result="${WORK_DIRECTORY}/frame_${state_dimension}.txt"

    if ! "${COMPILER_COMMAND}" "${BASE_FLAGS[@]}" "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" \
        "-DWATERMARK_STATE_DIMENSION=${state_dimension}" \
        "-DWATERMARK_INPUT_DIMENSION=${INPUT_DIMENSION}" \
        -fstack-usage -c "${INSTRUMENT_SOURCE}" -o "${object}" >"${diagnostics}" 2>&1
    then
        printf 'compile-failed\tcompile-failed\tcompile-failed\n' >"${result}"
        return
    fi

    local deepest_line
    deepest_line="$(grep -v 'stack_watermark\.cpp' "${report}" | sort -t$'\t' -k2 -nr | head -1)"

    local deepest_bytes deepest_function
    deepest_bytes="$(printf '%s' "${deepest_line}" | cut -f2)"
    deepest_function="$(printf '%s' "${deepest_line}" | cut -f1 | sed 's/.*:[0-9]*:[0-9]*://' |
        sed 's/(.*//' | sed 's/.* //')"

    local named_bytes
    named_bytes="$(grep "${PUBLISHED_FRAME_FUNCTION}" "${report}" | cut -f2 | sort -nr | head -1)"

    printf '%s\t%s\t%s\n' "${deepest_bytes:-none}" "${deepest_function:-none}" "${named_bytes:-none}" \
        >"${result}"
}

running=0
for state_dimension in "${STATE_DIMENSIONS[@]}"; do
    measure_dimension "${state_dimension}" &
    running=$((running + 1))
    if [ "${running}" -ge "${PARALLEL_JOBS}" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

if [ "${WITH_FRAMES}" -eq 1 ]; then
    running=0
    for state_dimension in "${STATE_DIMENSIONS[@]}"; do
        measure_frames "${state_dimension}" &
        running=$((running + 1))
        if [ "${running}" -ge "${PARALLEL_JOBS}" ]; then
            wait -n
            running=$((running - 1))
        fi
    done
    wait
fi

read_field()
{
    printf '%s' "$1" | tr ' ' '\n' | grep "^$2=" | cut -d= -f2-
}

printf 'Whole-chain runtime stack watermark of the discrete Riccati solve, swept\n'
printf 'over the state dimension, against the figures the real-time safety matrix\n'
printf 'publishes. The comparison is exact and carries no tolerance.\n\n'

printf 'Provenance\n'
printf '  compiler          %s\n' "${COMPILER_VERSION}"
printf '  flags             %s\n' "${BASE_FLAGS[*]}"
printf '  architecture flag none (compiler driver default)\n'
printf '  host              %s %s\n' "${HOST_SYSTEM}" "${HOST_ARCHITECTURE}"
printf '  scalar            double\n'
printf '  thread stack      %s MiB, set through the thread attribute\n' "${THREAD_STACK_MIB}"
printf '  corpus            discrete damped chain, input dimension %s\n' "${INPUT_DIMENSION}"
printf '\n'

FAILURE_COUNT=0
NONZERO_FLOOR_COUNT=0

printf '%-4s %-4s %-12s %-8s %-12s %-9s %-6s\n' \
    NX NU measured floor published delta match
for state_dimension in "${STATE_DIMENSIONS[@]}"; do
    line="$(cat "${WORK_DIRECTORY}/watermark_${state_dimension}.txt")"
    measured="$(read_field "${line}" watermark_bytes)"
    floor="$(read_field "${line}" floor_bytes)"
    published="${PUBLISHED_PEAK[${state_dimension}]}"

    if [ -z "${measured}" ]; then
        printf '%-4s %-4s %-12s %-8s %-12s %-9s %-6s\n' \
            "${state_dimension}" "${INPUT_DIMENSION}" "${line}" "--" "${published}" "--" "no"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
        continue
    fi

    delta=$((measured - published))
    if [ "${measured}" -eq "${published}" ]; then
        match="yes"
    else
        match="no"
        FAILURE_COUNT=$((FAILURE_COUNT + 1))
    fi
    if [ "${floor}" -ne 0 ]; then
        NONZERO_FLOOR_COUNT=$((NONZERO_FLOOR_COUNT + 1))
    fi

    printf '%-4s %-4s %-12s %-8s %-12s %-9s %-6s\n' \
        "${state_dimension}" "${INPUT_DIMENSION}" "${measured}" "${floor}" \
        "${published}" "${delta}" "${match}"
done
printf '\n'

if [ "${WITH_FRAMES}" -eq 1 ]; then
    printf 'Per-function frames. The report attributes nothing to callees, so these are\n'
    printf 'a LOWER BOUND and do not substitute for the watermark above. The published\n'
    printf 'column carries the acceptance check function, which is not in general the\n'
    printf 'deepest frame in the chain.\n\n'
    printf '%-4s %-10s %-46s %-14s %-12s %-6s\n' \
        NX deepest owner acceptance published match
    for state_dimension in "${STATE_DIMENSIONS[@]}"; do
        IFS=$'\t' read -r deepest_bytes deepest_function named_bytes \
            <"${WORK_DIRECTORY}/frame_${state_dimension}.txt"
        published="${PUBLISHED_FRAME[${state_dimension}]}"
        if [ "${named_bytes}" = "${published}" ]; then
            match="yes"
        else
            match="no"
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        printf '%-4s %-10s %-46s %-14s %-12s %-6s\n' \
            "${state_dimension}" "${deepest_bytes}" "${deepest_function}" \
            "${named_bytes}" "${published}" "${match}"
    done
    printf '\n'
fi

# One machine-readable record per dimension, so a consumer never has to parse the
# table above, and so the harness floor travels on the same line as the figure it
# qualifies. A watermark quoted without its floor is not a measurement.
for state_dimension in "${STATE_DIMENSIONS[@]}"; do
    line="$(cat "${WORK_DIRECTORY}/watermark_${state_dimension}.txt")"
    measured="$(read_field "${line}" watermark_bytes)"
    frame_bytes="--"
    if [ "${WITH_FRAMES}" -eq 1 ]; then
        frame_bytes="$(cut -f3 "${WORK_DIRECTORY}/frame_${state_dimension}.txt")"
    fi
    printf 'record row=riccati-discrete nx=%s nu=%s measured_peak_bytes=%s acceptance_frame_bytes=%s harness_floor_bytes=%s\n' \
        "${state_dimension}" "${INPUT_DIMENSION}" "${measured:-none}" "${frame_bytes}" \
        "$(read_field "${line}" floor_bytes)"
done
printf '\n'

if [ "${NONZERO_FLOOR_COUNT}" -ne 0 ]; then
    printf 'HARNESS FLOOR IS NONZERO on %s configuration(s). The figures above are not\n' \
        "${NONZERO_FLOOR_COUNT}"
    printf 'attributable to the call under measurement and nothing may be derived from\n'
    printf 'them until the floor reads zero.\n'
fi

if [ "${FAILURE_COUNT}" -eq 0 ]; then
    printf 'POSITIVE CONTROL PASSED: every figure equals the published one exactly.\n'
    exit 0
fi

printf 'POSITIVE CONTROL FAILED on %s cell(s).\n\n' "${FAILURE_COUNT}"
printf 'The instrument is NOT adjusted to make these agree, and no figure derived from\n'
printf 'it may be published until the discrepancy is attributed. The candidate\n'
printf 'mechanisms are the corpus, the optimization level, an architecture flag from\n'
printf 'the environment, a harness frame shifting the origin, and a change in the\n'
printf 'library since the figures were taken. The provenance block above pins the\n'
printf 'middle three; the corpus and the harness shape are the two the written record\n'
printf 'of the original measurement does not preserve.\n'
exit 1
