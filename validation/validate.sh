#!/usr/bin/env bash
# validate.sh -- Configure, build and run the whole cross-validation suite.
#
# Usage: ./validate.sh [build_dir]
#   build_dir  : CMake build directory (default: ./build)
#
# Every case is a test, and this script's exit status is the test runner's:
# a case whose executable was never built is a failure, not a printed line. To
# run a subset, invoke the runner directly with its own selection flags, e.g.
# `ctest --test-dir build -R dare_solution`.
#
# Each case produces an analysis/ subdirectory with CSVs, plots, and a report.
#
# Environment:
#   CTRLPP_VALIDATE_JOBS       build parallelism (default: 2)
#   CTRLPP_VALIDATE_GENERATOR  CMake generator (default: Unix Makefiles)
#   CTRLPP_VALIDATE_OCTAVE     interpreter to invoke (default: octave)
#
# Prerequisites:
#   - the interpreter named above, with the control, signal, splines, and
#     quaternion packages

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${1:-${SCRIPT_DIR}/build}"

JOBS="${CTRLPP_VALIDATE_JOBS:-2}"
GENERATOR="${CTRLPP_VALIDATE_GENERATOR:-Unix Makefiles}"
OCTAVE="${CTRLPP_VALIDATE_OCTAVE:-octave}"

echo "=== Building C++ validation cases ==="
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -G "$GENERATOR" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build "$BUILD_DIR" -j"$JOBS" 2>&1 | tail -5
echo ""

echo "=== Running validation cases ==="
set +e
ctest --test-dir "$BUILD_DIR" --no-tests=error --output-on-failure -T Test
STATUS=$?
set -e
echo ""

echo "=== Generating summary reports ==="
"$OCTAVE" --no-gui "$SCRIPT_DIR/generate_summary.m" "$SCRIPT_DIR" 2>/dev/null
echo ""

exit "$STATUS"
