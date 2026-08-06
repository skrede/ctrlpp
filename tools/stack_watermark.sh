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
# Usage: tools/stack_watermark.sh [--rows riccati|estimators|controller]
#                                 [--control | --diagonal | --held-dimension
#                                  | --interior | --ceilings]
#                                 [--frames] [--input-dimension N]
#                                 [--eigen PATH] [--backend-tree PATH]
#                                 [--gap-bytes N] [--paint-mib N]
#                                 [--stack-mib N] [--journal PATH] [--jobs N]
#
#   --rows riccati       the discrete Riccati rows, checked and unchecked
#   --rows estimators    the five estimator rows
#   --rows controller    the predictive controller row, whose build needs the
#                        optional nonlinear-programming backend
#   --control            the Riccati positive control (implies --rows riccati)
#   --diagonal           equal dimensions on every swept axis, with the
#                        per-grid-point build cost timed on a quiet station.
#                        The controller row takes its frame measurement here
#                        without being asked, because its published diagonal
#                        carries the frame column.
#   --held-dimension     each swept axis in turn with the others held, per row
#   --interior           the interior fill, over every row at once: a stride of
#                        four on every caller-chosen axis, from the first
#                        instantiable value to the instantiation limit. It is a
#                        build campaign rather than a build job, so it JOURNALS
#                        each finished point and skips what the journal already
#                        holds when it is re-run.
#   --ceilings           the two limits per row: the last dimension that
#                        compiles at all, walked to from the fill's top rung
#                        with the refusing diagnostic classified, and the last
#                        that fits each task stack, read off the fill's journal
#   --frames             also take the per-function frame measurement, which
#                        costs a second compile per grid point
#   --input-dimension N  the corpus input dimension (default 1); on the
#                        controller row the input dimension is a swept axis and
#                        this sets the value the other two sweeps hold it at
#   --eigen PATH         the linear-algebra include root (default the system one)
#   --backend-tree PATH  a cmake tree carrying the fetched nonlinear-programming
#                        backend. Configured FROM EMPTY at the backend's pinned
#                        revision when it is absent, so the figures belong to
#                        that revision rather than to a local working checkout,
#                        and the RESOLVED commit is read out of the fetched
#                        source rather than out of the pin that selected it.
#   --gap-bytes N        the harness gap, i.e. the measurement's resolution floor
#   --paint-mib N        the painted window in mebibytes. It is a CEILING on
#                        what the instrument can report, and the deepest chains
#                        in the interior of the grid reach the committed
#                        default, so the two stages that go there raise it.
#   --stack-mib N        the measurement thread's stack in mebibytes, which the
#                        instrument requires to be sixteen times the window
#   --journal PATH       where the interior fill records finished points, and
#                        where the limit stage reads the watermarks it turns
#                        into supported maxima
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
ROW_NMPC_STATIC=8

ESTIMATOR_ROWS=("${ROW_KALMAN}" "${ROW_EKF}" "${ROW_UKF}" "${ROW_MANIFOLD_UKF}" "${ROW_MEKF}")
RICCATI_ROWS=("${ROW_RICCATI_DISCRETE}" "${ROW_RICCATI_DISCRETE_UNCHECKED}")
CONTROLLER_ROWS=("${ROW_NMPC_STATIC}")

declare -A ROW_NAME=(
    ["${ROW_RICCATI_DISCRETE}"]=riccati-discrete
    ["${ROW_RICCATI_DISCRETE_UNCHECKED}"]=riccati-discrete-unchecked
    ["${ROW_KALMAN}"]=kalman-filter
    ["${ROW_EKF}"]=ekf
    ["${ROW_UKF}"]=ukf
    ["${ROW_MANIFOLD_UKF}"]=manifold-ukf
    ["${ROW_MEKF}"]=mekf
    ["${ROW_NMPC_STATIC}"]=nmpc-static
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
    ["${ROW_NMPC_STATIC}"]=state-dimension
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

# The predictive controller row's three CHOSEN dimensions, swept one at a time
# with the other two held, plus a diagonal. The values held are the
# configuration the shipped allocation proof pins, so the row's two pieces of
# evidence meet at one point.
#
# The ladders double while doubling is affordable and then bracket the
# instantiation ceiling, which on this row arrives at modest chosen dimensions:
# the dominant fixed-size object is the `NV x NV` Hessian, so the ceiling is a
# statement about the DERIVED decision dimension and the ladders approach it
# from three different directions. Each ladder's last two rungs sit either side
# of it, so the ceiling is bracketed rather than asserted.
# The horizon a row that has no horizon is compiled at. Only the controller row
# reads the selector, so the value is arbitrary and is named rather than spelled
# as a bare number at every call site.
HORIZON_NOT_READ=1

CONTROLLER_HELD_STATE=2
CONTROLLER_HELD_INPUT=1
CONTROLLER_HELD_HORIZON=5

CONTROLLER_STATE_LADDER=(1 2 4 8 16 20 21)
CONTROLLER_INPUT_LADDER=(1 2 4 8 16 23 24)
CONTROLLER_HORIZON_LADDER=(1 2 4 8 16 32 42 43)
CONTROLLER_DIAGONAL_LADDER=(1 2 3 4 5 6 7 8)

# The interior fill's axis: the smallest instantiable value, then every fourth
# value to the instantiation limit.
#
# ##########################################################################
# # THE STRIDE IS FOUR AND THE FILL THEREFORE DESCRIBES THE INTERIOR ON A   #
# # THIRTY-THREE-VALUE AXIS RATHER THAN EVERYWHERE.                         #
# #                                                                         #
# # Between two adjacent values on this axis nothing is measured, and no    #
# # claim is made about the three integers that lie there. What the stride  #
# # buys over sweeping each axis at one held value is the INTERACTION:      #
# # every value of one axis is measured against every value of the other,   #
# # so the shape of the dependence is described rather than only located.   #
# # What it does not buy is resolution between strides.                     #
# ##########################################################################
#
# The limit is 128 on every estimator axis, which is where a fixed-size
# `double` object of that dimension squared meets the linear-algebra library's
# fixed-size allocation limit exactly. The multiplicative filter's axis is the
# BIAS dimension and its error state is three larger, so its last stride value
# is 124 rather than 128 and its own limit is found by the limit stage rather
# than assumed here.
INTERIOR_AXIS=(1 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88
    92 96 100 104 108 112 116 120 124 128)
INTERIOR_BIAS_AXIS=(4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84
    88 92 96 100 104 108 112 116 120 124)

# The decision dimension the predictive controller row cannot exceed. Its
# ladders bracketed it from three directions in the plan that measured the row,
# and the interior fill enumerates the reachable triples under it rather than
# gridding a rectangle most of which does not instantiate.
CONTROLLER_DECISION_LIMIT=128

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
# What the stage calls itself on its summary line. It is a separate name from
# the selector because the summary line sits in the same stream as the per-row
# records, and a stage name that reads as one of the quantities those records
# carry would be counted as a row by anything scanning them.
STAGE_LABEL=""
ROW_SET="riccati"
GAP_BYTES=1024
EIGEN_ROOT="/usr/include/eigen3"
BACKEND_TREE=""
BACKEND_INCLUDE=""
BACKEND_REVISION="not-resolved"

# The committed defaults of the instrument's two window selectors, repeated
# here so the provenance block prints what was actually compiled with. They are
# CHECKED against the instrument below rather than trusted, because a default
# that drifts in one file and not the other publishes a window nothing was
# measured under.
PAINT_MIB=4
STACK_MIB=64
PAINT_MIB_EXPLICIT=0

# The window the two stages that reach the interior raise it to. The deepest
# chains there are within a factor of two of the committed 4 MiB window, and a
# chain that reaches the bottom of the window reports the window rather than
# itself. Eight times the deepest predicted chain is the headroom, and the
# saturation flag on every record is what makes the choice checkable rather
# than merely careful.
INTERIOR_PAINT_MIB=32
INTERIOR_STACK_MIB=512

JOURNAL=""
JOURNAL_PATH=""
DEFAULT_JOURNAL_NAME="stack_watermark_interior.records"

# How far above a value known to compile the limit stage is willing to walk
# before reporting that it did not find the boundary. A walk that runs off the
# end says so rather than reporting the last value it tried as a limit.
LIMIT_WALK_STEPS=16

# The task stacks the supported-maximum tables are published over.
TASK_STACK_LADDER=(4096 8192 16384 32768 49152 65536)

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
        --interior)
            STAGE="interior"
            shift
            ;;
        --ceilings)
            STAGE="ceilings"
            STAGE_LABEL="instantiation-and-stack-limits"
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
        --backend-tree)
            BACKEND_TREE="$2"
            shift 2
            ;;
        --gap-bytes)
            GAP_BYTES="$2"
            shift 2
            ;;
        --paint-mib)
            PAINT_MIB="$2"
            PAINT_MIB_EXPLICIT=1
            shift 2
            ;;
        --stack-mib)
            STACK_MIB="$2"
            PAINT_MIB_EXPLICIT=1
            shift 2
            ;;
        --journal)
            JOURNAL_PATH="$2"
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
    riccati|estimators|controller) ;;
    *)
        printf 'unrecognized row set: %s\n' "${ROW_SET}" >&2
        exit 2
        ;;
esac

if [ "${STAGE}" = "control" ] && [ "${ROW_SET}" != "riccati" ]; then
    printf 'the positive control is a Riccati measurement and has no form on other rows\n' >&2
    exit 2
fi

# The interior fill and the limit stage span every published row at once,
# including the one whose build needs the optional backend, so the row selector
# does not apply to them and is not silently honored.
SPANS_EVERY_ROW=0
if [ "${STAGE}" = "interior" ] || [ "${STAGE}" = "ceilings" ]; then
    SPANS_EVERY_ROW=1
    ROW_SET="every"
    if [ "${PAINT_MIB_EXPLICIT}" -eq 0 ]; then
        PAINT_MIB="${INTERIOR_PAINT_MIB}"
        STACK_MIB="${INTERIOR_STACK_MIB}"
    fi
fi

if ! command -v "${COMPILER_COMMAND}" >/dev/null 2>&1; then
    printf '%s is not available; nothing was measured\n' "${COMPILER_COMMAND}" >&2
    exit 1
fi

INCLUDE_FLAGS=(-I "${REPOSITORY_ROOT}/lib/ctrlpp/include" -isystem "${EIGEN_ROOT}")

# The optional nonlinear-programming backend, resolved for the controller row
# alone.
#
# ##########################################################################
# # THE PIN AND THE RESOLVED COMMIT ARE TWO DIFFERENT FACTS, AND ONLY THE   #
# # SECOND IS WHAT WAS MEASURED.                                           #
# #                                                                        #
# # The tree is configured FROM EMPTY so the backend is fetched at the      #
# # revision the project pins, rather than read out of whatever local       #
# # working checkout happens to be on the station carrying whatever         #
# # uncommitted work. The commit below is then read out of the FETCHED      #
# # source. Frames are an output of inlining that backend's templates, so a #
# # figure published without the revision beside it is a figure a reader    #
# # cannot re-derive.                                                      #
# ##########################################################################
resolve_backend()
{
    if [ -z "${BACKEND_TREE}" ]; then
        BACKEND_TREE="${REPOSITORY_ROOT}/build/stack_watermark_backend"
    fi

    local source_directory="${BACKEND_TREE}/_deps/argmin-src"
    if [ ! -d "${source_directory}" ]; then
        printf 'configuring %s from empty to fetch the backend at its pinned revision\n' "${BACKEND_TREE}" >&2
        rm -rf "${BACKEND_TREE}"
        if ! cmake -S "${REPOSITORY_ROOT}" -B "${BACKEND_TREE}" -G "Unix Makefiles" \
            -DCTRLPP_BUILD_ARGMIN=ON -DCTRLPP_CMAKE_FETCH_DEPS=ON \
            >"${WORK_DIRECTORY}/backend_configure.log" 2>&1
        then
            printf 'the backend could not be fetched; nothing was measured\n' >&2
            exit 1
        fi
    fi

    BACKEND_INCLUDE="${source_directory}/lib/argmin/include"
    if [ ! -d "${BACKEND_INCLUDE}" ]; then
        printf 'the fetched backend carries no include root at %s\n' "${BACKEND_INCLUDE}" >&2
        exit 1
    fi

    BACKEND_REVISION="$(git -C "${source_directory}" rev-parse HEAD 2>/dev/null)"
    BACKEND_REVISION="${BACKEND_REVISION:-unresolved}"
    INCLUDE_FLAGS+=(-isystem "${BACKEND_INCLUDE}")
}

# Named for the sweep rather than for the instrument, so that the directory does
# not collide with a single-configuration binary a developer builds by hand at
# the obvious path.
WORK_DIRECTORY="${TMPDIR:-/tmp}/stack_watermark_sweep"
rm -rf "${WORK_DIRECTORY}"
mkdir -p "${WORK_DIRECTORY}"

BACKEND_RESOLVED=0
if [ "${ROW_SET}" = "controller" ] || [ "${SPANS_EVERY_ROW}" -eq 1 ]; then
    resolve_backend
    BACKEND_RESOLVED=1
fi

COMPILER_VERSION="$("${COMPILER_COMMAND}" --version 2>/dev/null | head -1)"
HOST_ARCHITECTURE="$(uname -m)"
HOST_SYSTEM="$(uname -s)"

read_macro()
{
    grep -oE "define $1 [0-9]+" "${EIGEN_ROOT}/Eigen/src/Core/util/Macros.h" 2>/dev/null |
        grep -oE '[0-9]+' | head -1
}
EIGEN_VERSION="$(read_macro EIGEN_WORLD_VERSION).$(read_macro EIGEN_MAJOR_VERSION).$(read_macro EIGEN_MINOR_VERSION)"

# The instrument's own committed defaults for the painted window and the
# measurement thread's stack, read out of the instrument so this file cannot
# drift away from what is compiled. A drift is fatal rather than reported: the
# provenance block would otherwise name a window nothing was measured under.
read_default_mib()
{
    grep -oE "define $1 [0-9]+" "${INSTRUMENT_SOURCE}" | grep -oE '[0-9]+' | head -1
}

INSTRUMENT_PAINT_MIB="$(read_default_mib WATERMARK_PAINTED_MIB)"
INSTRUMENT_STACK_MIB="$(read_default_mib WATERMARK_THREAD_STACK_MIB)"

if [ "${INSTRUMENT_PAINT_MIB:-none}" != "4" ] || [ "${INSTRUMENT_STACK_MIB:-none}" != "64" ]; then
    printf 'the instrument default window (%s / %s MiB) is not the one this driver names\n' \
        "${INSTRUMENT_PAINT_MIB:-none}" "${INSTRUMENT_STACK_MIB:-none}" >&2
    exit 2
fi

if [ "$((STACK_MIB))" -lt "$((16 * PAINT_MIB))" ]; then
    printf 'a %s MiB window needs at least a %s MiB thread stack; %s was asked for\n' \
        "${PAINT_MIB}" "$((16 * PAINT_MIB))" "${STACK_MIB}" >&2
    exit 2
fi

THREAD_STACK_MIB="${STACK_MIB}"

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

# A grid point's own output path. Every dimension the instrument reads is in the
# name, because the controller row sweeps three of them and two points that
# differ only in the third would otherwise share a binary.
point_stem()
{
    printf '%s/%s_%s_%s_%s_%s' "${WORK_DIRECTORY}" "$1" "$2" "$3" "$4" "$5"
}

# The two DERIVED dimensions of the controller row, computed here from the
# relation as it is written in the published prose. The instrument computes them
# independently from its own template parameters and prints them, and the two
# are compared at every point: a disagreement means the relation is misstated in
# one of the two places, which is exactly the error a reader would inherit.
derived_decision_dimension()
{
    printf '%s' "$(( ($3 + 1) * $1 + $3 * $2 ))"
}

derived_constraint_dimension()
{
    printf '%s' "$(( $1 * ($3 + 1) ))"
}

# One grid point: compile and run in the same step, timing the compile and
# recording the station's condition on both sides of it. A wall-clock compile
# time taken beside another build measures station load rather than build cost,
# so the condition travels with the number instead of being asserted afterwards.
measure_point()
{
    local row="$1" state_dimension="$2" measurement_dimension="$3"
    local input_dimension="$4" horizon="$5" tag="$6"
    local stem
    stem="$(point_stem "${row}" "${state_dimension}" "${measurement_dimension}" \
        "${input_dimension}" "${horizon}")"
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
        "-DWATERMARK_INPUT_DIMENSION=${input_dimension}" \
        "-DWATERMARK_HORIZON=${horizon}" \
        "-DWATERMARK_HARNESS_GAP_BYTES=${GAP_BYTES}" \
        "-DWATERMARK_PAINTED_MIB=${PAINT_MIB}" \
        "-DWATERMARK_THREAD_STACK_MIB=${STACK_MIB}" \
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
    local input_dimension="$4" horizon="$5"
    local stem
    stem="$(point_stem "${row}" "${state_dimension}" "${measurement_dimension}" \
        "${input_dimension}" "${horizon}")"
    local object="${stem}.o"
    local report="${stem}.su"
    local diagnostics="${stem}.frame_diagnostics"
    local result="${stem}.frame"

    if ! "${COMPILER_COMMAND}" "${BASE_FLAGS[@]}" "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" \
        "-DWATERMARK_ROW=${row}" \
        "-DWATERMARK_STATE_DIMENSION=${state_dimension}" \
        "-DWATERMARK_MEASUREMENT_DIMENSION=${measurement_dimension}" \
        "-DWATERMARK_INPUT_DIMENSION=${input_dimension}" \
        "-DWATERMARK_HORIZON=${horizon}" \
        "-DWATERMARK_HARNESS_GAP_BYTES=${GAP_BYTES}" \
        "-DWATERMARK_PAINTED_MIB=${PAINT_MIB}" \
        "-DWATERMARK_THREAD_STACK_MIB=${STACK_MIB}" \
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
    printf '  painted window    %s MiB, i.e. the depth above which a chain reports the window\n' "${PAINT_MIB}"
    printf '                    instead of itself; every record carries whether it did\n'
    printf '  linear algebra    Eigen %s at %s\n' "${EIGEN_VERSION}" "${EIGEN_ROOT}"
    printf '  corpus            damped chain (vector-state rows) and constant body rate on SO(3)\n'
    printf '                    (attitude rows), input dimension %s\n' "${INPUT_DIMENSION}"
    printf '  concurrency       compiles run at -j%s\n' "${PARALLEL_JOBS}"
    if [ "${BACKEND_RESOLVED}" -eq 1 ]; then
        printf '  backend           argmin at RESOLVED commit %s\n' "${BACKEND_REVISION}"
        printf '                    fetched from empty at the project pin into %s\n' "${BACKEND_TREE}"
        printf '                    figures are NOT comparable across revisions of it\n'
    fi
    printf '\n'
}

NONZERO_FLOOR_COUNT=0
FAILURE_COUNT=0
BUSY_TIMING_COUNT=0
RELATION_MISMATCH_COUNT=0
SATURATED_COUNT=0
NOT_INSTANTIABLE_COUNT=0
MEASURED_COUNT=0

# A timed compile is qualified by the station's condition on BOTH sides of it. A
# nonzero count on either side means the number beside it measures station load
# rather than build cost.
#
# ##########################################################################
# # THE GATE APPLIES AT -j1 AND AT NO OTHER PARALLELISM, BECAUSE ABOVE IT   #
# # THE COUNTER SEES THIS DRIVER'S OWN SIBLING COMPILES.                    #
# #                                                                         #
# # At -j1 a per-point time IS a build cost and the station has to be quiet #
# # for it to be one. Above -j1 the per-point time is a campaign wall time  #
# # under the driver's own concurrency and is NOT comparable to a -j1       #
# # figure; the counter then reports the campaign's own siblings and cannot #
# # separate them from a foreign build. The count is still recorded on      #
# # every line, and the station is censused on both sides of the whole      #
# # stage, so the quietness condition is stated rather than asserted.       #
# ##########################################################################
count_busy_sides()
{
    [ "${PARALLEL_JOBS}" = "1" ] || return 0
    if [ "$1" != "--" ] && [ "$1" != "0" ]; then
        BUSY_TIMING_COUNT=$((BUSY_TIMING_COUNT + 1))
    fi
    if [ "$2" != "--" ] && [ "$2" != "0" ]; then
        BUSY_TIMING_COUNT=$((BUSY_TIMING_COUNT + 1))
    fi
}

# Every record goes to the standard output, and to the journal as well when one
# is in use.
#
# ##########################################################################
# # THE JOURNAL IS WHAT MAKES A CAMPAIGN SURVIVE AN INTERRUPTION.           #
# #                                                                         #
# # A fill measured in thousands of points is a build campaign rather than  #
# # a build job. One that must restart from zero after a failure part-way   #
# # is a campaign that does not finish, so each point is appended as it     #
# # completes, together with an index line naming the point, and a re-run   #
# # measures only what the index does not already hold.                     #
# ##########################################################################
emit_line()
{
    printf '%s\n' "$1"
    if [ -n "${JOURNAL}" ]; then
        printf '%s\n' "$1" >>"${JOURNAL}"
    fi
}

journal_key()
{
    printf '%s %s %s %s %s' "$1" "$2" "$3" "$4" "$5"
}

journal_holds()
{
    [ -n "${JOURNAL}" ] || return 1
    [ -f "${JOURNAL}.index" ] || return 1
    grep -qxF "$(journal_key "$@")" "${JOURNAL}.index"
}

journal_mark()
{
    [ -n "${JOURNAL}" ] || return 0
    printf '%s\n' "$(journal_key "$@")" >>"${JOURNAL}.index"
}

# A finished point's binary is deleted as soon as its record exists. The fill
# compiles thousands of them and keeping every one would fill the work
# directory without any of them ever being read again; the record is the
# artifact and the binary is not.
discard_point()
{
    local stem
    stem="$(point_stem "$1" "$2" "$3" "$4" "$5")"
    rm -f "${stem}.bin" "${stem}.o" "${stem}.su"
}

# One machine-readable record per grid point. The harness floor travels on the
# same line as the figure it qualifies, and so does the dimension held fixed: a
# watermark quoted without its floor is not a measurement, and one quoted without
# its held dimension is a figure a reader will take for a general one.
emit_record()
{
    local row="$1" state_dimension="$2" measurement_dimension="$3"
    local input_dimension="$4" horizon="$5" held="$6"
    local stem
    stem="$(point_stem "${row}" "${state_dimension}" "${measurement_dimension}" \
        "${input_dimension}" "${horizon}")"
    local line measured floor reported_nx reported_ny reported_nu
    line="$(cat "${stem}.result")"

    # The controller row's derived dimensions travel on the record whether the
    # point instantiated or not: the decision dimension at the point where the
    # instantiation stops compiling is precisely what a ceiling is a statement
    # about, and it is computed from the chosen dimensions rather than read back
    # from a binary that does not exist.
    local derived=""
    local expected_nv="" expected_maxm=""
    if [ "${row}" = "${ROW_NMPC_STATIC}" ]; then
        expected_nv="$(derived_decision_dimension "${state_dimension}" "${input_dimension}" "${horizon}")"
        expected_maxm="$(derived_constraint_dimension "${state_dimension}" "${input_dimension}" "${horizon}")"
        derived=" nh=${horizon} nv=${expected_nv} maxm=${expected_maxm}"
    fi

    # The first solve's whole-chain peak travels beside the steady-state one on
    # the row that has both, because the two are peaks of DIFFERENT chains and
    # the deeper of them is what a task stack has to cover.
    local first_call=""
    if [ "${row}" = "${ROW_NMPC_STATIC}" ] && [ -f "${stem}.result" ]; then
        local reported_first reported_construction
        reported_first="$(read_field "$(cat "${stem}.result")" first_call_bytes)"
        reported_construction="$(read_field "$(cat "${stem}.result")" construction_bytes)"
        if [ -n "${reported_first}" ]; then
            first_call=" first_call_watermark_bytes=${reported_first} construction_watermark_bytes=${reported_construction}"
        fi
    fi

    local seconds="--" busy_before="--" busy_after="--" load_before="--"
    if [ -f "${stem}.timing" ]; then
        IFS=$'\t' read -r seconds busy_before busy_after load_before _ <"${stem}.timing"
    fi

    if [ "${line}" = "compile-failed" ] || [ "${line}" = "run-failed" ]; then
        emit_line "$(printf '%s nx=%s ny=%s nu=%s%s held=%s not-instantiable=%s compile_seconds=%s jobs=-j%s loadavg_before=%s quiet_before=%s quiet_after=%s' \
            "${ROW_NAME[${row}]}" "${state_dimension}" "${measurement_dimension}" \
            "${input_dimension}" "${derived}" "${held}" "${line}" \
            "${seconds}" "${PARALLEL_JOBS}" "${load_before}" "${busy_before}" "${busy_after}")"
        NOT_INSTANTIABLE_COUNT=$((NOT_INSTANTIABLE_COUNT + 1))
        count_busy_sides "${busy_before}" "${busy_after}"
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
    MEASURED_COUNT=$((MEASURED_COUNT + 1))

    # A chain that reached the bottom of the painted window reports the window
    # rather than itself, and the difference is invisible in the figure. The
    # instrument decides this and the driver counts it, because a saturated
    # figure is a lower bound and no supported maximum may be derived from one.
    local saturated
    saturated="$(read_field "${line}" saturated)"
    if [ "${saturated}" = "yes" ]; then
        SATURATED_COUNT=$((SATURATED_COUNT + 1))
    fi

    # The instrument computed the same two derived dimensions from its own
    # template parameters. Comparing them here is what keeps the relation the
    # document publishes and the relation the measurement was taken at from
    # drifting apart silently.
    if [ "${row}" = "${ROW_NMPC_STATIC}" ]; then
        local reported_nv reported_maxm
        reported_nv="$(read_field "${line}" nv)"
        reported_maxm="$(read_field "${line}" maxm)"
        if [ "${reported_nv}" != "${expected_nv}" ] || [ "${reported_maxm}" != "${expected_maxm}" ]; then
            printf 'THE STATED RELATION DISAGREES WITH THE INSTANTIATION at nx=%s nu=%s nh=%s: stated (%s, %s), instantiated (%s, %s)\n' \
                "${state_dimension}" "${input_dimension}" "${horizon}" \
                "${expected_nv}" "${expected_maxm}" "${reported_nv}" "${reported_maxm}"
            RELATION_MISMATCH_COUNT=$((RELATION_MISMATCH_COUNT + 1))
        fi
    fi

    local frame_bytes="--" frame_function="--"
    if [ "${WITH_FRAMES}" -eq 1 ] && [ -f "${stem}.frame" ]; then
        IFS=$'\t' read -r frame_bytes frame_function _ <"${stem}.frame"
    fi

    emit_line "$(printf '%s nx=%s ny=%s nu=%s%s held=%s watermark_bytes=%s%s floor_bytes=%s saturated=%s deepest_frame_bytes=%s deepest_frame_function=%s compile_seconds=%s jobs=-j%s loadavg_before=%s quiet_before=%s quiet_after=%s' \
        "${ROW_NAME[${row}]}" "${reported_nx}" "${reported_ny}" "${reported_nu}" "${derived}" "${held}" \
        "${measured}" "${first_call}" "${floor}" "${saturated}" "${frame_bytes}" "${frame_function}" \
        "${seconds}" "${PARALLEL_JOBS}" "${load_before}" "${busy_before}" "${busy_after}")"

    count_busy_sides "${busy_before}" "${busy_after}"
}

run_grid()
{
    local -n points="$1"
    local running=0
    local point row state_dimension measurement_dimension input_dimension horizon
    for point in "${points[@]}"; do
        IFS=' ' read -r row state_dimension measurement_dimension input_dimension horizon _ <<<"${point}"
        measure_point "${row}" "${state_dimension}" "${measurement_dimension}" \
            "${input_dimension}" "${horizon}" "${point}" &
        running=$((running + 1))
        if [ "${running}" -ge "${PARALLEL_JOBS}" ]; then
            wait -n
            running=$((running - 1))
        fi
    done
    wait

    if [ "${WITH_FRAMES}" -eq 1 ]; then
        for point in "${points[@]}"; do
            IFS=' ' read -r row state_dimension measurement_dimension input_dimension horizon _ <<<"${point}"
            measure_frames "${row}" "${state_dimension}" "${measurement_dimension}" \
                "${input_dimension}" "${horizon}"
        done
    fi
}

STAGE_START="$(date +%s.%N)"
STAGE_START_BUSY="$(station_busy_count)"

if [ "${STAGE}" = "control" ]; then
    WITH_FRAMES=1
    POINTS=()
    for row in "${RICCATI_ROWS[@]}"; do
        for state_dimension in "${CONTROL_DIMENSIONS[@]}"; do
            POINTS+=("${row} ${state_dimension} 1 ${INPUT_DIMENSION} ${HORIZON_NOT_READ}")
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
            stem="$(point_stem "${row}" "${state_dimension}" 1 "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}")"
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
        stem="$(point_stem "${ROW_RICCATI_DISCRETE}" "${state_dimension}" 1 "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}")"
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
            emit_record "${row}" "${state_dimension}" 1 "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}" \
                "input-dimension-${INPUT_DIMENSION}"
        done
    done
    printf '\n'
fi

if [ "${STAGE}" = "diagonal" ] && [ "${ROW_SET}" = "controller" ]; then
    # The published diagonal for this row carries the frame column, so the frame
    # compile is taken here without being asked for rather than left to a flag a
    # reproducer could forget.
    WITH_FRAMES=1
    POINTS=()
    for row in "${CONTROLLER_ROWS[@]}"; do
        for rung in "${CONTROLLER_DIAGONAL_LADDER[@]}"; do
            POINTS+=("${row} ${rung} 0 ${rung} ${rung}")
        done
    done

    printf 'The diagonal: the state, input and horizon dimensions equal, one line per\n'
    printf 'rung. Every line carries the three CHOSEN dimensions, the two DERIVED\n'
    printf 'dimensions they induce, the whole-chain runtime watermark, the harness floor\n'
    printf 'that qualifies it, the deepest single frame with the function that owns it,\n'
    printf 'and the compile that produced the point with the station condition it was\n'
    printf 'timed under.\n\n'
    printf 'The derived dimensions are NOT gridded over and are not a caller choice:\n'
    printf '  NV   = (NH + 1) * NX + NH * NU\n'
    printf '  MaxM = NX * (NH + 1)\n'
    printf 'Most pairs of them correspond to no reachable configuration, so a grid over\n'
    printf 'them would describe a surface a caller cannot reach.\n\n'
    print_provenance
    run_grid POINTS

    for row in "${CONTROLLER_ROWS[@]}"; do
        for rung in "${CONTROLLER_DIAGONAL_LADDER[@]}"; do
            emit_record "${row}" "${rung}" 0 "${rung}" "${rung}" "nothing-diagonal"
        done
    done
    printf '\n'
fi

if [ "${STAGE}" = "held" ] && [ "${ROW_SET}" = "controller" ]; then
    POINTS=()
    for row in "${CONTROLLER_ROWS[@]}"; do
        for rung in "${CONTROLLER_STATE_LADDER[@]}"; do
            POINTS+=("${row} ${rung} 0 ${CONTROLLER_HELD_INPUT} ${CONTROLLER_HELD_HORIZON}")
        done
        for rung in "${CONTROLLER_INPUT_LADDER[@]}"; do
            POINTS+=("${row} ${CONTROLLER_HELD_STATE} 0 ${rung} ${CONTROLLER_HELD_HORIZON}")
        done
        for rung in "${CONTROLLER_HORIZON_LADDER[@]}"; do
            POINTS+=("${row} ${CONTROLLER_HELD_STATE} 0 ${CONTROLLER_HELD_INPUT} ${rung}")
        done
    done

    printf 'The held-dimension sweeps: each of the three CHOSEN dimensions swept in turn\n'
    printf 'with the other two held. Every line names what it holds and the value it\n'
    printf 'holds it at, and carries the two derived dimensions the point induces. No\n'
    printf 'interior point is measured: every point below lies on one of the three\n'
    printf 'held-dimension lines.\n\n'
    print_provenance
    run_grid POINTS

    for row in "${CONTROLLER_ROWS[@]}"; do
        for rung in "${CONTROLLER_STATE_LADDER[@]}"; do
            emit_record "${row}" "${rung}" 0 "${CONTROLLER_HELD_INPUT}" "${CONTROLLER_HELD_HORIZON}" \
                "input-dimension-at-${CONTROLLER_HELD_INPUT}-horizon-at-${CONTROLLER_HELD_HORIZON}"
        done
        for rung in "${CONTROLLER_INPUT_LADDER[@]}"; do
            emit_record "${row}" "${CONTROLLER_HELD_STATE}" 0 "${rung}" "${CONTROLLER_HELD_HORIZON}" \
                "state-dimension-at-${CONTROLLER_HELD_STATE}-horizon-at-${CONTROLLER_HELD_HORIZON}"
        done
        for rung in "${CONTROLLER_HORIZON_LADDER[@]}"; do
            emit_record "${row}" "${CONTROLLER_HELD_STATE}" 0 "${CONTROLLER_HELD_INPUT}" "${rung}" \
                "state-dimension-at-${CONTROLLER_HELD_STATE}-input-dimension-at-${CONTROLLER_HELD_INPUT}"
        done
    done
    printf '\n'
fi

if [ "${STAGE}" = "diagonal" ] && [ "${ROW_SET}" != "controller" ]; then
    POINTS=()
    for row in "${ESTIMATOR_ROWS[@]}"; do
        for rung in "${ESTIMATOR_LADDER[@]}"; do
            POINTS+=("${row} ${rung} ${rung} ${INPUT_DIMENSION} ${HORIZON_NOT_READ}")
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
            emit_record "${row}" "${rung}" "${rung}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}" "${held}"
        done
    done
    printf '\n'
fi

if [ "${STAGE}" = "held" ] && [ "${ROW_SET}" != "controller" ]; then
    POINTS=()
    for row in "${ESTIMATOR_ROWS[@]}"; do
        for rung in "${ESTIMATOR_LADDER[@]}"; do
            POINTS+=("${row} ${HELD_RUNG} ${rung} ${INPUT_DIMENSION} ${HORIZON_NOT_READ}")
        done
        # The manifold filter's state is a rotation, so there is no state
        # dimension for a caller to sweep and none is swept. Emitting one would
        # publish five readings of one configuration as if they were a trend.
        if [ "${row}" != "${ROW_MANIFOLD_UKF}" ]; then
            for rung in "${ESTIMATOR_LADDER[@]}"; do
                POINTS+=("${row} ${rung} ${HELD_RUNG} ${INPUT_DIMENSION} ${HORIZON_NOT_READ}")
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
            emit_record "${row}" "${HELD_RUNG}" "${rung}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}" \
                "${ROW_STATE_AXIS[${row}]}-at-${held_at}"
        done
        if [ "${row}" != "${ROW_MANIFOLD_UKF}" ]; then
            for rung in "${ESTIMATOR_LADDER[@]}"; do
                emit_record "${row}" "${rung}" "${HELD_RUNG}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}" \
                    "measurement-dimension-at-${HELD_RUNG}"
            done
        else
            printf '%s no-second-sweep held=%s-at-3\n' \
                "${ROW_NAME[${row}]}" "${ROW_STATE_AXIS[${row}]}"
        fi
    done
    printf '\n'
fi

SKIPPED_COUNT=0
PENDING=()

# One chunk of the campaign: every pending point is compiled and run
# concurrently, then each is recorded, journaled and discarded. The barrier
# between chunks costs the difference between the slowest and the mean point in
# the chunk, which is a few percent when a chunk holds one row's points, and it
# buys a bounded work directory and a journal that is never written from two
# processes at once.
drain_pending()
{
    local point row nx ny nu nh held
    for point in "${PENDING[@]}"; do
        IFS=' ' read -r row nx ny nu nh held <<<"${point}"
        measure_point "${row}" "${nx}" "${ny}" "${nu}" "${nh}" "${point}" &
    done
    wait

    for point in "${PENDING[@]}"; do
        IFS=' ' read -r row nx ny nu nh held <<<"${point}"
        emit_record "${row}" "${nx}" "${ny}" "${nu}" "${nh}" "${held}"
        journal_mark "${row}" "${nx}" "${ny}" "${nu}" "${nh}"
        discard_point "${row}" "${nx}" "${ny}" "${nu}" "${nh}"
    done
    PENDING=()
}

run_campaign()
{
    local -n campaign="$1"
    local total="${#campaign[@]}"
    local done_points=0
    local point row nx ny nu nh held
    for point in "${campaign[@]}"; do
        IFS=' ' read -r row nx ny nu nh held <<<"${point}"
        done_points=$((done_points + 1))
        if journal_holds "${row}" "${nx}" "${ny}" "${nu}" "${nh}"; then
            SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
            continue
        fi
        PENDING+=("${point}")
        if [ "${#PENDING[@]}" -ge "${PARALLEL_JOBS}" ]; then
            drain_pending
            printf 'progress %s/%s elapsed=%s\n' "${done_points}" "${total}" \
                "$(elapsed_seconds "${STAGE_START}" "$(date +%s.%N)")" >&2
        fi
    done
    if [ "${#PENDING[@]}" -gt 0 ]; then
        drain_pending
    fi
}

if [ "${STAGE}" = "interior" ]; then
    JOURNAL="${JOURNAL_PATH:-${REPOSITORY_ROOT}/build/${DEFAULT_JOURNAL_NAME}}"
    mkdir -p "$(dirname "${JOURNAL}")"

    POINTS=()
    for row in "${ROW_KALMAN}" "${ROW_EKF}" "${ROW_UKF}"; do
        for state_dimension in "${INTERIOR_AXIS[@]}"; do
            for measurement_dimension in "${INTERIOR_AXIS[@]}"; do
                POINTS+=("${row} ${state_dimension} ${measurement_dimension} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} nothing-interior-fill")
            done
        done
    done

    # One axis, because the rotation state is not a caller's to choose.
    for measurement_dimension in "${INTERIOR_AXIS[@]}"; do
        POINTS+=("${ROW_MANIFOLD_UKF} 3 ${measurement_dimension} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} ${ROW_STATE_AXIS[${ROW_MANIFOLD_UKF}]}-at-3")
    done

    # The bias axis starts at four rather than one: below three the propagation
    # has no leading bias block to subtract and the filter does not instantiate.
    # Three is legal and is not on the stride, which the published tables say.
    for state_dimension in "${INTERIOR_BIAS_AXIS[@]}"; do
        for measurement_dimension in "${INTERIOR_AXIS[@]}"; do
            POINTS+=("${ROW_MEKF} ${state_dimension} ${measurement_dimension} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} nothing-interior-fill")
        done
    done

    # The predictive controller row is enumerated over the REACHABLE triples
    # rather than over a rectangle: its decision dimension is derived, and most
    # of the rectangle induces one past the instantiation limit.
    for state_dimension in "${INTERIOR_AXIS[@]}"; do
        for input_dimension in "${INTERIOR_AXIS[@]}"; do
            for horizon in "${INTERIOR_AXIS[@]}"; do
                decision="$(derived_decision_dimension "${state_dimension}" "${input_dimension}" "${horizon}")"
                if [ "${decision}" -le "${CONTROLLER_DECISION_LIMIT}" ]; then
                    POINTS+=("${ROW_NMPC_STATIC} ${state_dimension} 0 ${input_dimension} ${horizon} nothing-interior-fill")
                fi
            done
        done
    done

    printf 'The interior fill: every caller-chosen axis of every row swept at a stride of\n'
    printf 'four, against every value of every other axis of that row. This describes the\n'
    printf 'shape of the interaction between the axes rather than only locating it, and it\n'
    printf 'resolves nothing between two adjacent values of a stride.\n\n'
    printf 'Points: %s. The run journals each finished point to\n' "${#POINTS[@]}"
    printf '%s and skips what that journal already holds, so an\n' "${JOURNAL}"
    printf 'interruption costs the chunk in flight rather than the campaign.\n\n'
    print_provenance
    run_campaign POINTS
    printf '\n'
fi

# The two limits, per row.
#
# ##########################################################################
# # THEY ANSWER DIFFERENT QUESTIONS AND NEITHER SUBSTITUTES FOR THE OTHER.  #
# #                                                                        #
# # The first is the last dimension that COMPILES: past it the             #
# # instantiation does not exist, and no task stack changes that. The      #
# # second is the last dimension whose whole-chain peak FITS a given task  #
# # stack: past it the instantiation exists and the task overflows. A      #
# # caller who does not fit a dimension needs to know which of the two it  #
# # is, because a bigger stack is the answer to one of them and nothing at #
# # all is the answer to the other.                                        #
# ##########################################################################
if [ "${STAGE}" = "ceilings" ]; then
    FILL_JOURNAL="${JOURNAL_PATH:-${REPOSITORY_ROOT}/build/${DEFAULT_JOURNAL_NAME}}"
    POINTS=()

    # A compile alone, with no link and no run: the question is whether the
    # instantiation exists, and a refusal is a diagnostic rather than a figure.
    # The refusal is CLASSIFIED, because an unrelated instantiation error at the
    # same dimension would look identical from the exit status and is not the
    # limit being reported.
    probe_instantiation()
    {
        local row="$1" state_dimension="$2" measurement_dimension="$3"
        local input_dimension="$4" horizon="$5"
        local log="${WORK_DIRECTORY}/probe_${row}_${state_dimension}_${measurement_dimension}_${input_dimension}_${horizon}.log"

        if "${COMPILER_COMMAND}" "${BASE_FLAGS[@]}" "${WARNING_FLAGS[@]}" "${INCLUDE_FLAGS[@]}" \
            "-DWATERMARK_ROW=${row}" \
            "-DWATERMARK_STATE_DIMENSION=${state_dimension}" \
            "-DWATERMARK_MEASUREMENT_DIMENSION=${measurement_dimension}" \
            "-DWATERMARK_INPUT_DIMENSION=${input_dimension}" \
            "-DWATERMARK_HORIZON=${horizon}" \
            "-DWATERMARK_HARNESS_GAP_BYTES=${GAP_BYTES}" \
            "-DWATERMARK_PAINTED_MIB=${PAINT_MIB}" \
            "-DWATERMARK_THREAD_STACK_MIB=${STACK_MIB}" \
            -fsyntax-only "${INSTRUMENT_SOURCE}" >"${log}" 2>&1
        then
            printf 'compiles'
        elif grep -q 'OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG' "${log}"; then
            printf 'refused-fixed-size-allocation-limit'
        else
            printf 'refused-unrelated-diagnostic'
        fi
    }

    # The tuple a value of the named axis stands for.
    point_for_axis()
    {
        local axis="$1" value="$2" other="$3"
        case "${axis}" in
            equal)
                printf '%s %s %s %s' "${value}" "${value}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}"
                ;;
            state)
                printf '%s %s %s %s' "${value}" "${other}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}"
                ;;
            measurement)
                printf '%s %s %s %s' "${other}" "${value}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}"
                ;;
            rotation-measurement)
                printf '3 %s %s %s' "${value}" "${INPUT_DIMENSION}" "${HORIZON_NOT_READ}"
                ;;
            controller-state)
                printf '%s 0 %s %s' "${value}" "${CONTROLLER_HELD_INPUT}" "${CONTROLLER_HELD_HORIZON}"
                ;;
            controller-input)
                printf '%s 0 %s %s' "${CONTROLLER_HELD_STATE}" "${value}" "${CONTROLLER_HELD_HORIZON}"
                ;;
            controller-horizon)
                printf '%s 0 %s %s' "${CONTROLLER_HELD_STATE}" "${CONTROLLER_HELD_INPUT}" "${value}"
                ;;
        esac
    }

    # Walk upward from a value the fill already compiled until the first
    # refusal. It is a WALK rather than a prediction: the starting value is one
    # the campaign above built and ran, and the boundary is the first value that
    # stops building. Two outcomes are reported per axis, one on each side of
    # the boundary, because a single refusal establishes that something failed
    # there rather than that the boundary is where it was expected.
    walk_axis()
    {
        local row="$1" axis="$2" start="$3" other="$4" held="$5"
        local value="${start}" verdict tuple
        local below="none" refused="none" refusal="none" steps=0

        tuple="$(point_for_axis "${axis}" "${value}" "${other}")"
        # shellcheck disable=SC2086
        verdict="$(probe_instantiation "${row}" ${tuple})"
        if [ "${verdict}" != "compiles" ]; then
            printf '%s axis=%s held=%s ceiling_kind=instantiation start=%s START_DOES_NOT_COMPILE=%s\n' \
                "${ROW_NAME[${row}]}" "${axis}" "${held}" "${start}" "${verdict}"
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
            return
        fi
        below="${value}"

        while [ "${steps}" -lt "${LIMIT_WALK_STEPS}" ]; do
            value=$((value + 1))
            steps=$((steps + 1))
            tuple="$(point_for_axis "${axis}" "${value}" "${other}")"
            # shellcheck disable=SC2086
            verdict="$(probe_instantiation "${row}" ${tuple})"
            if [ "${verdict}" = "compiles" ]; then
                below="${value}"
                continue
            fi
            refused="${value}"
            refusal="${verdict}"
            break
        done

        if [ "${refused}" = "none" ]; then
            printf '%s axis=%s held=%s ceiling_kind=instantiation last_compiling=%s BOUNDARY_NOT_REACHED_WITHIN=%s\n' \
                "${ROW_NAME[${row}]}" "${axis}" "${held}" "${below}" "${LIMIT_WALK_STEPS}"
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
            return
        fi

        printf '%s axis=%s held=%s ceiling_kind=instantiation last_compiling=%s first_refused=%s refusal=%s method=walked\n' \
            "${ROW_NAME[${row}]}" "${axis}" "${held}" "${below}" "${refused}" "${refusal}"

        if [ "${refusal}" != "refused-fixed-size-allocation-limit" ]; then
            printf 'THE REFUSAL AT %s ON %s IS NOT THE FIXED-SIZE ALLOCATION LIMIT, so the value\n' \
                "${refused}" "${ROW_NAME[${row}]}"
            printf 'above is where something else stopped the build and is not the limit this\n'
            printf 'stage reports.\n'
            FAILURE_COUNT=$((FAILURE_COUNT + 1))
        fi
    }

    printf 'The two limits per row. The first is the last dimension that compiles at all,\n'
    printf 'walked to from the interior fill top rung one integer at a time, with the\n'
    printf 'refusing diagnostic classified: an unrelated instantiation error at the same\n'
    printf 'dimension would look identical from the exit status alone. The second is the\n'
    printf 'last dimension whose whole-chain peak fits a given task stack, read off the\n'
    printf 'fill journal rather than off any per-function frame.\n\n'
    print_provenance

    printf 'The first limit, walked:\n\n'
    for row in "${ROW_KALMAN}" "${ROW_EKF}" "${ROW_UKF}"; do
        walk_axis "${row}" equal 128 128 nothing
        walk_axis "${row}" state 128 "${HELD_RUNG}" "measurement-dimension-at-${HELD_RUNG}"
        walk_axis "${row}" measurement 128 "${HELD_RUNG}" "state-dimension-at-${HELD_RUNG}"
    done
    walk_axis "${ROW_MANIFOLD_UKF}" rotation-measurement 128 3 "rotation-state-at-3"
    walk_axis "${ROW_MEKF}" equal 124 124 nothing
    walk_axis "${ROW_MEKF}" state 124 "${HELD_RUNG}" "measurement-dimension-at-${HELD_RUNG}"
    walk_axis "${ROW_MEKF}" measurement 128 "${HELD_RUNG}" "bias-dimension-at-${HELD_RUNG}"
    walk_axis "${ROW_NMPC_STATIC}" controller-state 20 0 \
        "input-dimension-at-${CONTROLLER_HELD_INPUT}-horizon-at-${CONTROLLER_HELD_HORIZON}"
    walk_axis "${ROW_NMPC_STATIC}" controller-input 23 0 \
        "state-dimension-at-${CONTROLLER_HELD_STATE}-horizon-at-${CONTROLLER_HELD_HORIZON}"
    walk_axis "${ROW_NMPC_STATIC}" controller-horizon 42 0 \
        "state-dimension-at-${CONTROLLER_HELD_STATE}-input-dimension-at-${CONTROLLER_HELD_INPUT}"
    printf '\n'

    # One measured point per row at the value the walk confirmed, so this stage
    # emits a watermark with its harness floor beside it rather than resting
    # entirely on compiles that were never run.
    POINTS=(
        "${ROW_KALMAN} ${HELD_RUNG} ${HELD_RUNG} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} nothing-limit-witness"
        "${ROW_EKF} ${HELD_RUNG} ${HELD_RUNG} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} nothing-limit-witness"
        "${ROW_UKF} ${HELD_RUNG} ${HELD_RUNG} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} nothing-limit-witness"
        "${ROW_MANIFOLD_UKF} 3 ${HELD_RUNG} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} rotation-state-at-3"
        "${ROW_MEKF} ${HELD_RUNG} ${HELD_RUNG} ${INPUT_DIMENSION} ${HORIZON_NOT_READ} nothing-limit-witness"
        "${ROW_NMPC_STATIC} ${CONTROLLER_HELD_STATE} 0 ${CONTROLLER_HELD_INPUT} ${CONTROLLER_HELD_HORIZON} nothing-limit-witness"
    )
    run_grid POINTS
    printf 'A measured point per row, so every line of this stage carries a harness floor:\n\n'
    for point in "${POINTS[@]}"; do
        IFS=' ' read -r row state_dimension measurement_dimension input_dimension horizon held <<<"${point}"
        emit_record "${row}" "${state_dimension}" "${measurement_dimension}" \
            "${input_dimension}" "${horizon}" "${held}"
    done
    printf '\n'

    printf 'The second limit, from the whole-chain watermarks of the fill:\n\n'
    if [ ! -f "${FILL_JOURNAL}" ]; then
        printf 'NO FILL JOURNAL AT %s, so no supported maximum is reported here.\n' "${FILL_JOURNAL}"
        printf 'The per-function frame is NOT substituted for it: that report attributes\n'
        printf 'nothing to callees and would overstate what fits a task stack.\n\n'
    else
        for stack_bytes in "${TASK_STACK_LADDER[@]}"; do
            awk -v limit="${stack_bytes}" '
                # The rotation-state row holds its state at three rather than at
                # the four the other rows hold theirs at, because three is what
                # it structurally is. Reading its swept line at four would find
                # nothing and report a zero that means "not measured" where the
                # other rows report one that means "does not fit".
                function held_state(row) { return row == "manifold-ukf" ? 3 : 4 }

                /watermark_bytes=/ && /saturated=no/ {
                    row = $1
                    nx = ny = nv = -1
                    peak = first_solve = -1
                    for(field = 2; field <= NF; ++field)
                    {
                        split($field, pair, "=")
                        if(pair[1] == "nx") nx = pair[2] + 0
                        else if(pair[1] == "ny") ny = pair[2] + 0
                        else if(pair[1] == "nv") nv = pair[2] + 0
                        else if(pair[1] == "watermark_bytes") peak = pair[2] + 0
                        else if(pair[1] == "first_call_watermark_bytes") first_solve = pair[2] + 0
                    }
                    if(peak < 0) next

                    # A caller that solves at all runs the first solve once, and
                    # on the row that has a first solve deeper than its steady
                    # state a maximum taken from the steady state alone would be
                    # an understatement. An understated stack figure overflows a
                    # task rather than returning a wrong answer.
                    if(first_solve > peak) peak = first_solve
                    seen[row] = 1

                    # The predictive controller row is not gridded over two
                    # named axes: its three chosen dimensions induce a decision
                    # dimension, and the peak is dominated by that. The honest
                    # statement over a fill is therefore the decision dimension
                    # below which EVERY measured configuration fits, reported
                    # beside the deepest single one that does.
                    if(row == "nmpc-static")
                    {
                        if(nv < 0) next
                        if(peak <= limit)
                        {
                            if(nv > deepest_fitting[row]) deepest_fitting[row] = nv
                        }
                        else if(refused[row] == 0 || nv < refused[row])
                            refused[row] = nv
                        next
                    }

                    if(peak > limit) next
                    if(nx == ny && nx > equal[row]) equal[row] = nx
                    if(nx == held_state(row) && ny > swept_measurement[row]) swept_measurement[row] = ny
                    if(ny == 4 && nx > swept_state[row]) swept_state[row] = nx
                }
                END {
                    for(row in seen)
                    {
                        if(row == "nmpc-static")
                            printf "%s task_stack_bytes=%d largest_fitting_decision_dimension=%d every_measured_configuration_fits_below=%d source=whole-chain-watermark\n",
                                row, limit, deepest_fitting[row], refused[row]
                        else if(row == "manifold-ukf")
                            printf "%s task_stack_bytes=%d largest_fitting_measurement_axis=%d rotation-state-is-not-a-caller-axis=yes source=whole-chain-watermark\n",
                                row, limit, swept_measurement[row]
                        else
                            printf "%s task_stack_bytes=%d largest_fitting_equal=%d largest_fitting_measurement_axis=%d largest_fitting_state_axis=%d source=whole-chain-watermark\n",
                                row, limit, equal[row], swept_measurement[row], swept_state[row]
                    }
                }
            ' "${FILL_JOURNAL}" | sort
        done
        printf '\n'
        printf 'A zero above means no value on that line of that row fits that task stack.\n'
        printf 'The measurement-axis column holds the state axis at four, or at the\n'
        printf 'structural three on the rotation row; the state-axis column holds the\n'
        printf 'measurement axis at four. Every figure is the whole-chain runtime\n'
        printf 'watermark; the per-function frame is a lower bound that excludes every\n'
        printf 'library frame beneath the named function and is never substituted for it.\n\n'
    fi
fi

STAGE_FINISH="$(date +%s.%N)"
printf 'stage=%s points=%s measured=%s not_instantiable=%s already_journaled=%s cumulative_wall_seconds=%s jobs=-j%s loadavg_at_end=%s quiet_at_stage_start=%s quiet_at_stage_end=%s\n\n' \
    "${STAGE_LABEL:-${STAGE}}" "${#POINTS[@]}" "${MEASURED_COUNT}" "${NOT_INSTANTIABLE_COUNT}" "${SKIPPED_COUNT}" \
    "$(elapsed_seconds "${STAGE_START}" "${STAGE_FINISH}")" \
    "${PARALLEL_JOBS}" "$(one_minute_load)" "${STAGE_START_BUSY}" "$(station_busy_count)"

if [ "${SATURATED_COUNT}" -ne 0 ]; then
    printf 'THE PAINTED WINDOW SATURATED on %s configuration(s). Those chains reached the\n' \
        "${SATURATED_COUNT}"
    printf 'bottom of the window, so the walk reported where it stopped looking rather than\n'
    printf 'where the chain stopped writing. Each such figure is a lower bound and no\n'
    printf 'supported maximum may be derived from one; re-run those points with a larger\n'
    printf 'window rather than publishing them.\n'
    exit 1
fi

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

if [ "${RELATION_MISMATCH_COUNT}" -ne 0 ]; then
    printf 'THE PUBLISHED RELATION AND THE INSTANTIATION DISAGREE on %s point(s). The\n' \
        "${RELATION_MISMATCH_COUNT}"
    printf 'derived dimensions above are computed twice, once from the relation as it is\n'
    printf 'written and once from the template parameters the measurement ran at, and a\n'
    printf 'reader would inherit whichever of the two is wrong.\n'
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
