#!/usr/bin/env bash
#
# Driver for the stack watermark instrument beside this file.
#
# It compiles and runs `stack_watermark.cpp` once per grid point, under the
# provenance every figure it prints is published with, and prints the whole-chain
# runtime watermark and the harness floor on the same line.
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
# Each grid point is compiled and run in the same step, to its own output path.
# No binary is reused across grid points. A table whose whole content is the
# difference between grid points says nothing at all if a stale binary answers
# for one of them, and it still looks complete.
#
# ##########################################################################
# # EVERY FIGURE CARRIES ITS PROVENANCE, AND THE INPUT DIMENSION AND THE   #
# # LINEAR-ALGEBRA RELEASE ARE PART OF IT.                                 #
# #                                                                        #
# # The input dimension moves the discrete Riccati watermark by more than a #
# # thousand bytes at two states, and the linear-algebra library's patch    #
# # release moves other measured quantities in this project by the same     #
# # magnitude as fused-multiply-add contraction does. Both are selectable   #
# # here, both are printed in the provenance block, and both travel on the  #
# # machine-readable record lines. A figure quoted without them cannot be   #
# # re-derived, which is how a published stack table stops being            #
# # reproducible.                                                          #
# ##########################################################################
#
# The comparison against the recorded control figures is EXACT and carries no
# tolerance. A harness adjusted until the figures agree is not a control, so a
# mismatch is printed as a mismatch and the driver exits nonzero. So does a
# nonzero harness floor, which means the figure beside it is not attributable to
# the call under measurement.
#
# Usage: tools/stack_watermark.sh [--rows riccati|estimators]
#                                 [--control | --diagonal | --held-dimension]
#                                 [--frames] [--input-dimension N]
#                                 [--eigen PATH] [--gap-bytes N] [--jobs N]
#
#   --rows riccati       the discrete Riccati rows, checked and unchecked
#   --rows estimators    the five estimator rows
#   --control            the Riccati positive control (implies --rows riccati)
#   --diagonal           equal state and measurement dimension, with the
#                        per-grid-point build cost timed on a quiet station
#   --held-dimension     measurement swept at fixed state, then state swept at
#                        fixed measurement, per row
#   --frames             also take the per-function frame measurement, which
#                        costs a second compile per grid point
#   --input-dimension N  the corpus input dimension (default 1)
#   --eigen PATH         the linear-algebra include root (default the system one)
#   --gap-bytes N        the harness gap, i.e. the measurement's resolution floor
#   --jobs N             concurrent compiles (default 1; the build-cost timings
#                        are only meaningful at 1)

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
# promotions, and the linear-algebra headers as a system include exactly as the
# library's build tree treats them, so a diagnostic reported below comes from the
# instrument's own source.
WARNING_FLAGS=(
    -fPIC -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
    -Wold-style-cast -Wcast-align -Woverloaded-virtual -Wnon-virtual-dtor
    -Wdouble-promotion -Wimplicit-fallthrough -Wformat=2
    -Werror=unused-result -Werror=unused-value
)

# The row selectors, which must agree with the instrument's own.
ROW_RICCATI_DISCRETE=1
ROW_RICCATI_DISCRETE_UNCHECKED=2
ROW_KALMAN=3
ROW_EKF=4
ROW_UKF=5
ROW_MANIFOLD_UKF=6
ROW_MEKF=7

ESTIMATOR_ROWS=("${ROW_KALMAN}" "${ROW_EKF}" "${ROW_UKF}" "${ROW_MANIFOLD_UKF}" "${ROW_MEKF}")
RICCATI_ROWS=("${ROW_RICCATI_DISCRETE}" "${ROW_RICCATI_DISCRETE_UNCHECKED}")

declare -A ROW_NAME=(
    ["${ROW_RICCATI_DISCRETE}"]=riccati-discrete
    ["${ROW_RICCATI_DISCRETE_UNCHECKED}"]=riccati-discrete-unchecked
    ["${ROW_KALMAN}"]=kalman-filter
    ["${ROW_EKF}"]=ekf
    ["${ROW_UKF}"]=ukf
    ["${ROW_MANIFOLD_UKF}"]=manifold-ukf
    ["${ROW_MEKF}"]=mekf
)

# What a row's first axis IS. It is not "the state dimension" in three of the
# five cases, and a held-dimension label that says it is would be the exact error
# the held-dimension requirement exists to prevent: the multiplicative filter's
# selector is the BIAS dimension and its error state is three larger, and the
# manifold filter has no such selector at all because its state is a rotation.
declare -A ROW_STATE_AXIS=(
    ["${ROW_RICCATI_DISCRETE}"]=state-dimension
    ["${ROW_RICCATI_DISCRETE_UNCHECKED}"]=state-dimension
    ["${ROW_KALMAN}"]=state-dimension
    ["${ROW_EKF}"]=state-dimension
    ["${ROW_UKF}"]=state-dimension
    ["${ROW_MANIFOLD_UKF}"]=rotation-state-dimension-not-caller-controlled
    ["${ROW_MEKF}"]=bias-dimension
)

# The state dimension the Riccati control is taken at, and the geometric ladder
# the estimator stages are taken over. The ladder doubles, so a quadratic cost
# quadruples between rungs and the growth is readable off four steps.
CONTROL_DIMENSIONS=(2 4 6 8)
ESTIMATOR_LADDER=(2 4 8 16 32)

# The rung every held-dimension sweep holds its other dimension at. It is the
# smallest rung at which every one of the five rows instantiates: the
# multiplicative filter's bias dimension is at least three, so the ladder's first
# rung does not exist for it.
HELD_RUNG=4

# The whole-chain runtime peaks and the acceptance-check frames the real-time
# safety matrix publishes for the discrete Riccati rows, keyed by state
# dimension, at input dimension 1 against the system linear-algebra release.
# These are the values THIS instrument measures on THIS tree; the control asks
# that an extension to the instrument leaves them alone.
declare -A CONTROL_PEAK_CHECKED=([2]=4552 [4]=11176 [6]=21912 [8]=41304)
declare -A CONTROL_PEAK_UNCHECKED=([2]=3368 [4]=8176 [6]=13560 [8]=25416)
declare -A CONTROL_FRAME=([2]=864 [4]=2400 [6]=6432 [8]=15056)

# The published frame column carries THIS function's frame, and it is not in
# general the deepest frame in the chain: at two states a Schur reordering
# routine has a larger one. The driver prints both rather than letting either
# stand for the other.
CONTROL_FRAME_FUNCTION="estimate_riccati_forward_error"

INPUT_DIMENSION=1
PARALLEL_JOBS=1
WITH_FRAMES=0
STAGE="control"
ROW_SET="riccati"
GAP_BYTES=1024
EIGEN_ROOT="/usr/include/eigen3"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --frames)
            WITH_FRAMES=1
            shift
            ;;
        --control)
            STAGE="control"
            ROW_SET="riccati"
            shift
            ;;
        --diagonal)
            STAGE="diagonal"
            shift
            ;;
        --held-dimension)
            STAGE="held"
            shift
            ;;
        --rows)
            ROW_SET="$2"
            shift 2
            ;;
        --input-dimension)
            INPUT_DIMENSION="$2"
            shift 2
            ;;
        --eigen)
            EIGEN_ROOT="$2"
            shift 2
            ;;
        --gap-bytes)
            GAP_BYTES="$2"
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

case "${ROW_SET}" in
    riccati|estimators) ;;
    *)
        printf 'unrecognized row set: %s\n' "${ROW_SET}" >&2
        exit 2
        ;;
esac

if [ "${STAGE}" = "control" ] && [ "${ROW_SET}" = "estimators" ]; then
    printf 'the positive control is a Riccati measurement and has no estimator form\n' >&2
    exit 2
fi

if ! command -v "${COMPILER_COMMAND}" >/dev/null 2>&1; then
    printf '%s is not available; nothing was measured\n' "${COMPILER_COMMAND}" >&2
    exit 1
fi

INCLUDE_FLAGS=(-I "${REPOSITORY_ROOT}/lib/ctrlpp/include" -isystem "${EIGEN_ROOT}")

# Named for the sweep rather than for the instrument, so that the directory does
# not collide with a single-configuration binary a developer builds by hand at
# the obvious path.
WORK_DIRECTORY="${TMPDIR:-/tmp}/stack_watermark_sweep"
rm -rf "${WORK_DIRECTORY}"
mkdir -p "${WORK_DIRECTORY}"

COMPILER_VERSION="$("${COMPILER_COMMAND}" --version 2>/dev/null | head -1)"
HOST_ARCHITECTURE="$(uname -m)"
HOST_SYSTEM="$(uname -s)"

read_macro()
{
    grep -oE "define $1 [0-9]+" "${EIGEN_ROOT}/Eigen/src/Core/util/Macros.h" 2>/dev/null |
        grep -oE '[0-9]+' | head -1
}
EIGEN_VERSION="$(read_macro EIGEN_WORLD_VERSION).$(read_macro EIGEN_MAJOR_VERSION).$(read_macro EIGEN_MINOR_VERSION)"

# The thread stack the instrument requests, read out of the instrument itself so
# the provenance block cannot drift away from what is actually measured.
THREAD_STACK_MIB="$(grep -oE 'thread_stack_bytes = std::size_t\{[0-9]+\}' "${INSTRUMENT_SOURCE}" |
    grep -oE '[0-9]+' | tail -1)"

# The station-quiet check. It prints the COUNT, and the count is what is
# examined: `pgrep` exits nonzero precisely when nothing matches, which here is
# the passing outcome, so the exit status is the wrong thing to test.
BUILD_PROCESS_PATTERN='cc1plus|cc1|cc|clang|clang\+\+|g\+\+|make|cmake|ld'

station_busy_count()
{
    # `pgrep -c` prints the count and exits NONZERO precisely when the count is
    # zero, which here is the passing outcome. The printed count is therefore
    # what is read; the exit status is discarded rather than tested.
    local count
    count="$(pgrep -c -x "${BUILD_PROCESS_PATTERN}" 2>/dev/null)"
    printf '%s' "${count:-0}"
}

one_minute_load()
{
    cut -d' ' -f1 /proc/loadavg 2>/dev/null || printf 'unavailable'
}

elapsed_seconds()
{
    awk -v start="$1" -v finish="$2" 'BEGIN { printf "%.2f", finish - start }'
}

# One grid point: compile and run in the same step, timing the compile and
# recording the station's condition on both sides of it. A wall-clock compile
# time taken beside another build measures station load rather than build cost,
# so the condition travels with the number instead of being asserted afterwards.
measure_point()
{
    local row="$1" state_dimension="$2" measurement_dimension="$3" tag="$4"
    local stem="${WORK_DIRECTORY}/${row}_${state_dimension}_${measurement_dimension}"
    local binary="${stem}.bin"
    local diagnostics="${stem}.diagnostics"
    local result="${stem}.result"
    local timing="${stem}.timing"

    local busy_before load_before start finish busy_after
    busy_before="$(station_busy_count)"
    load_before="$(one_minute_load)"
    start="$(date +%s.%N)"

    if "${COMPILER_COMMAND}" "${BASE_FLAGS[@]}" "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" \
        "-DWATERMARK_ROW=${row}" \
        "-DWATERMARK_STATE_DIMENSION=${state_dimension}" \
        "-DWATERMARK_MEASUREMENT_DIMENSION=${measurement_dimension}" \
        "-DWATERMARK_INPUT_DIMENSION=${INPUT_DIMENSION}" \
        "-DWATERMARK_HARNESS_GAP_BYTES=${GAP_BYTES}" \
        "${INSTRUMENT_SOURCE}" -o "${binary}" >"${diagnostics}" 2>&1
    then
        finish="$(date +%s.%N)"
        if ! "${binary}" >"${result}" 2>>"${diagnostics}"; then
            printf 'run-failed\n' >"${result}"
        fi
    else
        finish="$(date +%s.%N)"
        printf 'compile-failed\n' >"${result}"
    fi

    busy_after="$(station_busy_count)"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$(elapsed_seconds "${start}" "${finish}")" "${busy_before}" "${busy_after}" \
        "${load_before}" "${tag}" >"${timing}"
}

# The name of the function a frame report line belongs to.
#
# The signature the report prints carries the return type, the template
# arguments and the parameter list, and the template arguments THEMSELVES
# contain parentheses -- a dimension spelled as a cast is the common case. So
# cutting at the first parenthesis lands inside the return type and leaves no
# name at all. The template groups are removed innermost-first until none
# remain, which is what makes the first remaining parenthesis the parameter
# list's, and the name is then the last token before it.
owning_function()
{
    local text previous
    text="$(printf '%s' "$1" | sed 's/.*:[0-9]*:[0-9]*://')"
    while :; do
        previous="${text}"
        text="$(printf '%s' "${text}" | sed 's/<[^<>]*>//g')"
        [ "${text}" = "${previous}" ] && break
    done
    text="$(printf '%s' "${text}" | sed 's/(.*//' | sed 's/ *$//' | sed 's/.* //')"
    printf '%s' "${text:-none}"
}

# The per-function frame report. The compiler writes one line per EMITTED
# function, so a function the optimizer inlined has no line of its own and its
# frame is merged into its caller's. That is why the instrument's own
# translation unit is excluded from the "deepest" search below: the chain is
# inlined into the instrument's call wrapper at the larger dimensions, and
# reporting the wrapper would attribute a library frame to the harness.
measure_frames()
{
    local row="$1" state_dimension="$2" measurement_dimension="$3"
    local stem="${WORK_DIRECTORY}/${row}_${state_dimension}_${measurement_dimension}"
    local object="${stem}.o"
    local report="${stem}.su"
    local diagnostics="${stem}.frame_diagnostics"
    local result="${stem}.frame"

    if ! "${COMPILER_COMMAND}" "${BASE_FLAGS[@]}" "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" \
        "-DWATERMARK_ROW=${row}" \
        "-DWATERMARK_STATE_DIMENSION=${state_dimension}" \
        "-DWATERMARK_MEASUREMENT_DIMENSION=${measurement_dimension}" \
        "-DWATERMARK_INPUT_DIMENSION=${INPUT_DIMENSION}" \
        "-DWATERMARK_HARNESS_GAP_BYTES=${GAP_BYTES}" \
        -fstack-usage -c "${INSTRUMENT_SOURCE}" -o "${object}" >"${diagnostics}" 2>&1
    then
        printf 'compile-failed\tcompile-failed\tcompile-failed\n' >"${result}"
        return
    fi

    local deepest_line
    deepest_line="$(grep -v 'stack_watermark\.cpp' "${report}" | sort -t$'\t' -k2 -nr | head -1)"

    local deepest_bytes deepest_function
    deepest_bytes="$(printf '%s' "${deepest_line}" | cut -f2)"
    deepest_function="$(owning_function "$(printf '%s' "${deepest_line}" | cut -f1)")"

    local named_bytes
    named_bytes="$(grep "${CONTROL_FRAME_FUNCTION}" "${report}" | cut -f2 | sort -nr | head -1)"

    printf '%s\t%s\t%s\n' "${deepest_bytes:-none}" "${deepest_function:-none}" "${named_bytes:-none}" \
        >"${result}"
}

read_field()
{
    printf '%s' "$1" | tr ' ' '\n' | grep "^$2=" | cut -d= -f2-
}

print_provenance()
{
    printf 'Provenance\n'
    printf '  compiler          %s\n' "${COMPILER_VERSION}"
    printf '  flags             %s\n' "${BASE_FLAGS[*]}"
    printf '  architecture flag none (compiler driver default)\n'
    printf '  host              %s %s\n' "${HOST_SYSTEM}" "${HOST_ARCHITECTURE}"
    printf '  scalar            double\n'
    printf '  thread stack      %s MiB, set through the thread attribute\n' "${THREAD_STACK_MIB}"
    printf '  harness gap       %s bytes, i.e. the resolution below which a chain reports zero\n' "${GAP_BYTES}"
    printf '  linear algebra    Eigen %s at %s\n' "${EIGEN_VERSION}" "${EIGEN_ROOT}"
    printf '  corpus            damped chain (vector-state rows) and constant body rate on SO(3)\n'
    printf '                    (attitude rows), input dimension %s\n' "${INPUT_DIMENSION}"
    printf '  concurrency       compiles run at -j%s\n' "${PARALLEL_JOBS}"
    printf '\n'
}

NONZERO_FLOOR_COUNT=0
FAILURE_COUNT=0
BUSY_TIMING_COUNT=0

# One machine-readable record per grid point. The harness floor travels on the
# same line as the figure it qualifies, and so does the dimension held fixed: a
# watermark quoted without its floor is not a measurement, and one quoted without
# its held dimension is a figure a reader will take for a general one.
emit_record()
{
    local row="$1" state_dimension="$2" measurement_dimension="$3" held="$4"
    local stem="${WORK_DIRECTORY}/${row}_${state_dimension}_${measurement_dimension}"
    local line measured floor reported_nx reported_ny reported_nu
    line="$(cat "${stem}.result")"

    if [ "${line}" = "compile-failed" ] || [ "${line}" = "run-failed" ]; then
        printf '%s nx=%s ny=%s held=%s not-instantiable=%s\n' \
            "${ROW_NAME[${row}]}" "${state_dimension}" "${measurement_dimension}" "${held}" "${line}"
        return
    fi

    measured="$(read_field "${line}" watermark_bytes)"
    floor="$(read_field "${line}" floor_bytes)"
    reported_nx="$(read_field "${line}" nx)"
    reported_ny="$(read_field "${line}" ny)"
    reported_nu="$(read_field "${line}" nu)"
    if [ "${floor}" -ne 0 ]; then
        NONZERO_FLOOR_COUNT=$((NONZERO_FLOOR_COUNT + 1))
    fi

    local frame_bytes="--" frame_function="--"
    if [ "${WITH_FRAMES}" -eq 1 ] && [ -f "${stem}.frame" ]; then
        IFS=$'\t' read -r frame_bytes frame_function _ <"${stem}.frame"
    fi

    local seconds="--" busy_before="--" busy_after="--" load_before="--"
    if [ -f "${stem}.timing" ]; then
        IFS=$'\t' read -r seconds busy_before busy_after load_before _ <"${stem}.timing"
    fi

    printf '%s nx=%s ny=%s nu=%s held=%s watermark_bytes=%s floor_bytes=%s deepest_frame_bytes=%s deepest_frame_function=%s compile_seconds=%s jobs=-j%s loadavg_before=%s quiet_before=%s quiet_after=%s\n' \
        "${ROW_NAME[${row}]}" "${reported_nx}" "${reported_ny}" "${reported_nu}" "${held}" \
        "${measured}" "${floor}" "${frame_bytes}" "${frame_function}" \
        "${seconds}" "${PARALLEL_JOBS}" "${load_before}" "${busy_before}" "${busy_after}"

    if [ "${busy_before}" != "--" ] && [ "${busy_before}" != "0" ]; then
        BUSY_TIMING_COUNT=$((BUSY_TIMING_COUNT + 1))
    fi
    if [ "${busy_after}" != "--" ] && [ "${busy_after}" != "0" ]; then
        BUSY_TIMING_COUNT=$((BUSY_TIMING_COUNT + 1))
    fi
}

run_grid()
{
    local -n points="$1"
    local running=0
    local point row state_dimension measurement_dimension
    for point in "${points[@]}"; do
        IFS=' ' read -r row state_dimension measurement_dimension _ <<<"${point}"
        measure_point "${row}" "${state_dimension}" "${measurement_dimension}" "${point}" &
        running=$((running + 1))
        if [ "${running}" -ge "${PARALLEL_JOBS}" ]; then
            wait -n
            running=$((running - 1))
        fi
    done
    wait

    if [ "${WITH_FRAMES}" -eq 1 ]; then
        for point in "${points[@]}"; do
            IFS=' ' read -r row state_dimension measurement_dimension _ <<<"${point}"
            measure_frames "${row}" "${state_dimension}" "${measurement_dimension}"
        done
    fi
}

STAGE_START="$(date +%s.%N)"

if [ "${STAGE}" = "control" ]; then
    WITH_FRAMES=1
    POINTS=()
    for row in "${RICCATI_ROWS[@]}"; do
        for state_dimension in "${CONTROL_DIMENSIONS[@]}"; do
            POINTS+=("${row} ${state_dimension} 1")
        done
    done

    printf 'Whole-chain runtime stack watermark of the discrete Riccati solve, swept\n'
    printf 'over the state dimension, against the recorded control figures. The\n'
    printf 'comparison is exact and carries no tolerance.\n\n'
    print_provenance
    run_grid POINTS

    printf '%-28s %-4s %-4s %-12s %-8s %-12s %-9s %-6s\n' \
        row NX NU measured floor recorded delta match
    for row in "${RICCATI_ROWS[@]}"; do
        for state_dimension in "${CONTROL_DIMENSIONS[@]}"; do
            stem="${WORK_DIRECTORY}/${row}_${state_dimension}_1"
            line="$(cat "${stem}.result")"
            measured="$(read_field "${line}" watermark_bytes)"
            reported_floor="$(read_field "${line}" floor_bytes)"
            if [ "${row}" = "${ROW_RICCATI_DISCRETE}" ]; then
                recorded="${CONTROL_PEAK_CHECKED[${state_dimension}]}"
            else
                recorded="${CONTROL_PEAK_UNCHECKED[${state_dimension}]}"
            fi

            if [ -z "${measured}" ]; then
                printf '%-28s %-4s %-4s %-12s %-8s %-12s %-9s %-6s\n' \
                    "${ROW_NAME[${row}]}" "${state_dimension}" "${INPUT_DIMENSION}" \
                    "${line}" "--" "${recorded}" "--" "no"
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
                continue
            fi

            delta=$((measured - recorded))
            if [ "${measured}" -eq "${recorded}" ]; then
                match="yes"
            else
                match="no"
                FAILURE_COUNT=$((FAILURE_COUNT + 1))
            fi
            if [ "${reported_floor}" -ne 0 ]; then
                NONZERO_FLOOR_COUNT=$((NONZERO_FLOOR_COUNT + 1))
            fi

            printf '%-28s %-4s %-4s %-12s %-8s %-12s %-9s %-6s\n' \
                "${ROW_NAME[${row}]}" "${state_dimension}" "${INPUT_DIMENSION}" \
                "${measured}" "${reported_floor}" "${recorded}" "${delta}" "${match}"
        done
    done
    printf '\n'

    printf 'Per-function frames. The report attributes nothing to callees, so these are\n'
    printf 'a LOWER BOUND and do not substitute for the watermark above. The recorded\n'
    printf 'column carries the acceptance check function, which is not in general the\n'
    printf 'deepest frame in the chain.\n\n'
    printf '%-4s %-10s %-46s %-14s %-12s %-6s\n' \
        NX deepest owner acceptance recorded match
    for state_dimension in "${CONTROL_DIMENSIONS[@]}"; do
        stem="${WORK_DIRECTORY}/${ROW_RICCATI_DISCRETE}_${state_dimension}_1"
        IFS=$'\t' read -r deepest_bytes deepest_function named_bytes <"${stem}.frame"
        recorded="${CONTROL_FRAME[${state_dimension}]}"
        if [ "${named_bytes}" = "${recorded}" ]; then
            match="yes"
        else
            match="no"
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
        printf '%-4s %-10s %-46s %-14s %-12s %-6s\n' \
            "${state_dimension}" "${deepest_bytes}" "${deepest_function}" \
            "${named_bytes}" "${recorded}" "${match}"
    done
    printf '\n'

    for row in "${RICCATI_ROWS[@]}"; do
        for state_dimension in "${CONTROL_DIMENSIONS[@]}"; do
            emit_record "${row}" "${state_dimension}" 1 "input-dimension-${INPUT_DIMENSION}"
        done
    done
    printf '\n'
fi

if [ "${STAGE}" = "diagonal" ]; then
    POINTS=()
    for row in "${ESTIMATOR_ROWS[@]}"; do
        for rung in "${ESTIMATOR_LADDER[@]}"; do
            POINTS+=("${row} ${rung} ${rung}")
        done
    done

    printf 'Stage one, the diagonal: equal state and measurement dimension, one line per\n'
    printf 'row per rung. Every line carries the whole-chain runtime watermark, the\n'
    printf 'harness floor that qualifies it, the deepest single frame with the function\n'
    printf 'that owns it, and the compile that produced the point with the station\n'
    printf 'condition it was timed under.\n\n'
    print_provenance
    run_grid POINTS

    for row in "${ESTIMATOR_ROWS[@]}"; do
        held="nothing-diagonal"
        if [ "${row}" = "${ROW_MANIFOLD_UKF}" ]; then
            held="${ROW_STATE_AXIS[${row}]}-at-3"
        fi
        for rung in "${ESTIMATOR_LADDER[@]}"; do
            emit_record "${row}" "${rung}" "${rung}" "${held}"
        done
    done
    printf '\n'
fi

if [ "${STAGE}" = "held" ]; then
    POINTS=()
    for row in "${ESTIMATOR_ROWS[@]}"; do
        for rung in "${ESTIMATOR_LADDER[@]}"; do
            POINTS+=("${row} ${HELD_RUNG} ${rung}")
        done
        # The manifold filter's state is a rotation, so there is no state
        # dimension for a caller to sweep and none is swept. Emitting one would
        # publish five readings of one configuration as if they were a trend.
        if [ "${row}" != "${ROW_MANIFOLD_UKF}" ]; then
            for rung in "${ESTIMATOR_LADDER[@]}"; do
                POINTS+=("${row} ${rung} ${HELD_RUNG}")
            done
        fi
    done

    printf 'Stage two, the held-dimension sweeps: the measurement dimension swept at a\n'
    printf 'fixed state dimension, then the state dimension swept at a fixed measurement\n'
    printf 'dimension. Every line names the dimension it holds and the value it holds it\n'
    printf 'at. No interior point is measured: every point below lies on one of the two\n'
    printf 'held-dimension lines, and the point where they cross is on the diagonal.\n\n'
    print_provenance
    run_grid POINTS

    for row in "${ESTIMATOR_ROWS[@]}"; do
        held_at="${HELD_RUNG}"
        if [ "${row}" = "${ROW_MANIFOLD_UKF}" ]; then
            held_at=3
        fi
        for rung in "${ESTIMATOR_LADDER[@]}"; do
            emit_record "${row}" "${HELD_RUNG}" "${rung}" "${ROW_STATE_AXIS[${row}]}-at-${held_at}"
        done
        if [ "${row}" != "${ROW_MANIFOLD_UKF}" ]; then
            for rung in "${ESTIMATOR_LADDER[@]}"; do
                emit_record "${row}" "${rung}" "${HELD_RUNG}" "measurement-dimension-at-${HELD_RUNG}"
            done
        else
            printf '%s no-second-sweep held=%s-at-3\n' \
                "${ROW_NAME[${row}]}" "${ROW_STATE_AXIS[${row}]}"
        fi
    done
    printf '\n'
fi

STAGE_FINISH="$(date +%s.%N)"
printf 'stage=%s points=%s cumulative_wall_seconds=%s jobs=-j%s loadavg_at_end=%s\n\n' \
    "${STAGE}" "${#POINTS[@]}" "$(elapsed_seconds "${STAGE_START}" "${STAGE_FINISH}")" \
    "${PARALLEL_JOBS}" "$(one_minute_load)"

if [ "${BUSY_TIMING_COUNT}" -ne 0 ]; then
    printf 'THE STATION WAS NOT QUIET on %s side(s) of a timed compile. Those timings\n' \
        "${BUSY_TIMING_COUNT}"
    printf 'measure station load rather than build cost and must be re-taken on a quiet\n'
    printf 'station rather than annotated after the fact.\n\n'
fi

if [ "${NONZERO_FLOOR_COUNT}" -ne 0 ]; then
    printf 'HARNESS FLOOR IS NONZERO on %s configuration(s). The figures above are not\n' \
        "${NONZERO_FLOOR_COUNT}"
    printf 'attributable to the call under measurement and nothing may be derived from\n'
    printf 'them until the floor reads zero.\n'
    exit 1
fi

if [ "${STAGE}" != "control" ]; then
    exit 0
fi

if [ "${FAILURE_COUNT}" -eq 0 ]; then
    printf 'POSITIVE CONTROL PASSED: every figure equals the recorded one exactly.\n'
    exit 0
fi

printf 'POSITIVE CONTROL FAILED on %s cell(s).\n\n' "${FAILURE_COUNT}"
printf 'The instrument is NOT adjusted to make these agree, and no figure derived from\n'
printf 'it may be published until the discrepancy is attributed. The candidate\n'
printf 'mechanisms are the corpus and its input dimension, the linear-algebra release,\n'
printf 'the optimization level, an architecture flag from the environment, a harness\n'
printf 'frame shifting the origin, and a change in the library since the figures were\n'
printf 'taken. The provenance block above pins every one of them.\n'
exit 1
