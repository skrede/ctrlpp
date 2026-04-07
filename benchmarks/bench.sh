#!/usr/bin/env bash
# bench.sh -- Configure, build, and run ctrlpp benchmarks.
#
# Usage: ./bench.sh [build_dir] [--internal-only] [--comparison lib1,lib2,...]
#   build_dir       : CMake build directory (default: ./build)
#   --internal-only : Skip comparison benchmarks entirely
#   --comparison   : Enable specific competitors (comma-separated)
#                     Valid: libmpc,osqp_eigen,hpipm,ct,ruckig,drake,argmin

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
CMAKE_OPTS=""
INTERNAL_ONLY=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --internal-only)
            INTERNAL_ONLY=true
            shift
            ;;
        --comparison)
            shift
            IFS=',' read -ra COMPETITORS <<< "${1:-}"
            for comp in "${COMPETITORS[@]}"; do
                case "$comp" in
                    libmpc)     CMAKE_OPTS+=" -DCTRLPP_BENCH_LIBMPC=ON" ;;
                    osqp_eigen) CMAKE_OPTS+=" -DCTRLPP_BENCH_OSQP_EIGEN=ON" ;;
                    hpipm)      CMAKE_OPTS+=" -DCTRLPP_BENCH_HPIPM=ON" ;;
                    ct)         CMAKE_OPTS+=" -DCTRLPP_BENCH_CT=ON" ;;
                    ruckig)     CMAKE_OPTS+=" -DCTRLPP_BENCH_RUCKIG=ON" ;;
                    drake)      CMAKE_OPTS+=" -DCTRLPP_BENCH_DRAKE=ON" ;;
                    argmin)     CMAKE_OPTS+=" -DCTRLPP_BENCH_ARGMIN=ON" ;;
                    *)          echo "Unknown competitor: $comp" >&2; exit 1 ;;
                esac
            done
            shift
            ;;
        *)
            BUILD_DIR="$1"
            shift
            ;;
    esac
done

# Configure
echo "=== Configuring benchmarks ==="
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release ${CMAKE_OPTS} 2>&1 | tail -5
echo ""

# Build
echo "=== Building benchmarks ==="
cmake --build "$BUILD_DIR" -j"$(nproc)" 2>&1 | tail -5
echo ""

# Run benchmarks
echo "=== Running benchmarks ==="

BENCH_COUNT=0

run_benchmarks_in() {
    local dir="$1"
    if [[ ! -d "$dir" ]]; then
        return
    fi
    while IFS= read -r -d '' binary; do
        local name
        name="$(basename "$binary")"
        echo "--- $name ---"
        "$binary"
        echo ""
        BENCH_COUNT=$((BENCH_COUNT + 1))
    done < <(find "$dir" -maxdepth 2 -name 'bench_*' -executable -type f -print0 | sort -z)
}

# Internal benchmarks
run_benchmarks_in "$BUILD_DIR/internal"

# Competitive benchmarks (unless --internal-only)
if [[ "$INTERNAL_ONLY" == false ]]; then
    run_benchmarks_in "$BUILD_DIR/comparison"
fi

# Collect CSV results
RESULTS_DIR="$SCRIPT_DIR/results"
mkdir -p "$RESULTS_DIR"

csv_count=0
while IFS= read -r -d '' csv; do
    mv "$csv" "$RESULTS_DIR/"
    csv_count=$((csv_count + 1))
done < <(find "$BUILD_DIR" -name '*.csv' -type f -print0 2>/dev/null)

echo "=== Summary ==="
echo "Benchmarks run: $BENCH_COUNT"
if [[ $csv_count -gt 0 ]]; then
    echo "CSV results:    $csv_count files in $RESULTS_DIR/"
else
    echo "CSV results:    none (no benchmark targets built yet)"
fi
