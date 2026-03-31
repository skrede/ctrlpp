#!/usr/bin/env bash
# validate.sh -- Build C++ cases, run both Octave and C++ for each case, compare outputs.
#
# Usage: ./validate.sh [build_dir] [case_name...]
#   build_dir  : CMake build directory (default: ./build)
#   case_name  : run only named cases (default: all cases in cases/)
#
# Each case produces an analysis/ subdirectory with CSVs, plots, and a report.
#
# Prerequisites:
#   - octave-cli with control package
#   - C++ cases built via CMake (this script builds them)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${1:-${SCRIPT_DIR}/build}"
shift 2>/dev/null || true

# Build C++ cases
echo "=== Building C++ validation cases ==="
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build "$BUILD_DIR" -j"$(nproc)" 2>&1 | tail -5
echo ""

# Collect cases
if [ $# -gt 0 ]; then
    CASES=("$@")
else
    CASES=()
    for d in "$SCRIPT_DIR"/cases/*/; do
        [ -d "$d" ] && CASES+=("$(basename "$d")")
    done
fi

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for case_name in "${CASES[@]}"; do
    case_dir="$SCRIPT_DIR/cases/$case_name"
    octave_script="$case_dir/${case_name}.m"
    cpp_binary="$BUILD_DIR/cases/$case_name/$case_name"

    printf -- "--- %-40s " "$case_name"

    if [ ! -f "$octave_script" ]; then
        echo "SKIP (no .m file)"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    if [ ! -x "$cpp_binary" ]; then
        echo "SKIP (no C++ binary)"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    # Create analysis output directory
    analysis_dir="$case_dir/analysis"
    mkdir -p "$analysis_dir"

    ref_csv="$analysis_dir/${case_name}_octave.csv"
    cand_csv="$analysis_dir/${case_name}_cpp.csv"

    # Read tolerances from case config if present
    atol="1e-10"
    rtol="1e-8"
    if [ -f "$case_dir/tolerance.cfg" ]; then
        source "$case_dir/tolerance.cfg"
    fi

    # Run Octave reference
    if ! octave --no-gui "$octave_script" > "$ref_csv" 2>/dev/null; then
        echo "FAIL (octave error)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # Run C++ candidate
    if ! "$cpp_binary" > "$cand_csv" 2>/dev/null; then
        echo "FAIL (C++ error)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    echo ""

    # Compare, generate plots and report
    if octave --no-gui "$SCRIPT_DIR/validate_compare.m" "$ref_csv" "$cand_csv" "$atol" "$rtol" "$analysis_dir" 2>/dev/null; then
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    echo ""
done

echo "=== Summary: $PASS_COUNT passed, $FAIL_COUNT failed, $SKIP_COUNT skipped ==="
echo ""

# Generate top-level summary reports
echo "=== Generating summary reports ==="
octave --no-gui "$SCRIPT_DIR/generate_summary.m" "$SCRIPT_DIR" 2>/dev/null
echo ""

[ "$FAIL_COUNT" -eq 0 ]
