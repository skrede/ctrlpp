#!/usr/bin/env bash
#
# Enumerates the promotion and narrowing diagnostics across the single-precision
# translation units, one line per unit-and-compiler pair plus a total.
#
# Two diagnostics are counted: -Wdouble-promotion, which is the softfloat-timing
# and determinism hazard the single-precision surface exists to keep out, and
# the narrowing into single precision, which catches a decimal constant that is
# in range but not exactly representable there. They are counted together because
# neither compiler reports the whole union on its own -- gcc reports the
# conversions and clang reports the promotions -- so a single compiler's zero
# says nothing about the other's.
#
# The narrowing has two spellings. gcc files it under float-conversion together
# with a float-to-integer conversion; clang keeps float-conversion for the
# integer case only and reports the narrowing as implicit-float-conversion, which
# -Wconversion below already enables. Counting only gcc's spelling would make
# every clang narrowing invisible while still printing a zero.
#
# The warning set below mirrors the non-MSVC branch of CTRLPP_WARNING_FLAGS, and
# the include incantation mirrors what the build actually passes: Eigen arrives
# through -isystem, the test framework through -I, and the unit-test units carry
# the default tree's -fno-exceptions -fno-rtti and the framework's matching
# define. A looser incantation would either flood the output with third-party
# diagnostics or hide ones the build sees.
#
# This is a developer tool, deliberately not registered as a test: it drives two
# compilers over the same sources the build already compiles once each, and its
# value is the enumeration rather than a pass or fail in a suite.
#
# Usage: scripts/float_tier_diagnostics.sh [build-directory]
#
# The build directory supplies the fetched Eigen and test-framework headers and
# defaults to the tree the default configure preset produces.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

build_dir="${1:-build/dev}"

ctrlpp_include="lib/ctrlpp/include"
eigen_include="${build_dir}/_deps/eigen3-src"
catch_include="${build_dir}/_deps/catch2-src/src"
catch_generated="${build_dir}/_deps/catch2-build/generated-includes"

for required in "${eigen_include}/Eigen/Dense" \
                "${catch_include}/catch2/catch_test_macros.hpp" \
                "${catch_generated}/catch2/catch_user_config.hpp"
do
    if [ ! -f "${required}" ]; then
        echo "ERROR: ${required} is missing." >&2
        echo "  Pass a configured build directory as the first argument; the default is build/dev." >&2
        echo "  Any tree configured with CTRLPP_BUILD_TESTS=ON and CTRLPP_CMAKE_FETCH_DEPS=ON will do." >&2
        exit 1
    fi
done

warning_flags=(
    -fPIC
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wsign-conversion
    -Wold-style-cast
    -Wcast-align
    -Woverloaded-virtual
    -Wnon-virtual-dtor
    -Wdouble-promotion
    -Wfloat-conversion
    -Wimplicit-fallthrough
    -Wformat=2
)

units=(
    tests/compile/embedded_core_float.cpp
    tests/unit/float_instantiation_test.cpp
    tests/unit/float_tier_test.cpp
)

compilers=(g++ clang++)

counted='\[-W(double-promotion|float-conversion|implicit-float-conversion)\]'

total=0
printf '%-46s %-10s %s\n' "translation unit" "compiler" "diagnostics"

for unit in "${units[@]}"; do
    unit_args=(-I "${ctrlpp_include}" -isystem "${eigen_include}")
    case "${unit}" in
        tests/unit/*)
            # The default tree compiles its unit tests throw-free and hands the
            # test framework the matching define; without it the framework's own
            # headers stop matching the flags and the counts drift from the build.
            unit_args+=(-I "${catch_include}" -I "${catch_generated}"
                        -DCATCH_CONFIG_DISABLE_EXCEPTIONS -fno-exceptions -fno-rtti)
            ;;
    esac

    for compiler in "${compilers[@]}"; do
        if ! command -v "${compiler}" >/dev/null 2>&1; then
            printf '%-46s %-10s %s\n' "${unit}" "${compiler}" "SKIPPED (not on PATH)"
            continue
        fi

        output="$("${compiler}" -fsyntax-only -std=c++20 \
            "${unit_args[@]}" "${warning_flags[@]}" "${unit}" 2>&1)"
        hits="$(printf '%s\n' "${output}" | grep -E -c -- "${counted}")"
        total=$((total + hits))

        printf '%-46s %-10s %d\n' "${unit}" "${compiler}" "${hits}"
        if [ "${hits}" -gt 0 ]; then
            printf '%s\n' "${output}" | grep -E -- "${counted}" | sed 's/^/    /'
        fi
    done
done

printf '%-46s %-10s %d\n' "TOTAL" "" "${total}"

if [ "${total}" -ne 0 ]; then
    exit 1
fi
