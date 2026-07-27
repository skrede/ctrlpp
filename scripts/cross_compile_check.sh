#!/usr/bin/env bash
#
# Local embedded-clean cross-compile check for the core ctrlpp surface.
#
# Runs three legs, each of which must exit 0:
#   Leg 1  Host build with exceptions and RTTI off (-fno-exceptions -fno-rtti
#          -DCTRLPP_NO_EXCEPTIONS). Authoritative test that no throw or stray
#          .value() leaks onto the embedded path; compiled for Scalar double
#          and float.
#   Leg 2  Host -std=c++20 build with exceptions on and CTRLPP_NO_EXCEPTIONS
#          undefined, exercising the ctrlpp::detail::expected fallback and the
#          CTRLPP_HAS_EXCEPTIONS-gated wrappers in their compiled-in state. The
#          full build-tests suite is the deep fallback gate since the library
#          floor is C++20; this leg keeps the script self-contained.
#   Leg 3  arm-none-eabi Cortex-M7 cross-compile with exceptions and RTTI off,
#          proving the header subset is clean against a bare-metal toolchain.
#          The toolchain's own hosted libstdc++ subset is the target, so the
#          standard-library search path and hosting mode are left at their
#          defaults (see the D-13 recipe).
#
# This is a plain, repeatable local script. It authors no CI configuration and
# no CMake toolchain file.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

ctrlpp_include="lib/ctrlpp/include"
witness_tu="tests/compile/embedded_core_float.cpp"

host_cxx="${CXX:-g++}"
arm_cxx="arm-none-eabi-g++"

# --- Eigen include resolution -------------------------------------------------
# Prefer the FetchContent checkout under build-tests, then a system install.

eigen_include=""
for candidate in \
    "build-tests/_deps/eigen3-src" \
    "build-tests/_deps/eigen-src" \
    "/usr/include/eigen3"
do
    if [ -f "${candidate}/Eigen/Dense" ]; then
        eigen_include="${candidate}"
        break
    fi
done

if [ -z "${eigen_include}" ]; then
    echo "ERROR: could not locate the Eigen headers." >&2
    echo "  Looked for Eigen/Dense under:" >&2
    echo "    build-tests/_deps/eigen3-src" >&2
    echo "    build-tests/_deps/eigen-src" >&2
    echo "    /usr/include/eigen3" >&2
    echo "  Configure the test build (which fetches Eigen) or install Eigen >= 3.4." >&2
    exit 1
fi

echo "Using Eigen include: ${eigen_include}"
echo "Using ctrlpp include: ${ctrlpp_include}"
echo

# --- Optional double-instantiation probe --------------------------------------
# The committed witness is the float translation unit. To honor the double and
# float coverage of the host no-exceptions leg we additionally compile a small
# generated probe that instantiates the same core umbrella surface at double.

# mktemp -d is portable across GNU and BSD/macOS (unlike GNU-only --suffix),
# and a probe.cpp inside it keeps the .cpp extension the compiler needs to
# infer the C++ language front end.
probe_dir="$(mktemp -d)"
double_probe="${probe_dir}/probe.cpp"
trap 'rm -rf "${probe_dir}"' EXIT
cat > "${double_probe}" <<'PROBE'
#include "ctrlpp/control.h"
#include "ctrlpp/estimation.h"
#include "ctrlpp/dsp.h"
#include "ctrlpp/trajectory.h"

namespace
{
using scalar = double;

int fold_all()
{
    int witness = 0;

    ctrlpp::pid_config<scalar, 1> pid_cfg{};
    ctrlpp::pid<scalar, 1> controller(pid_cfg);
    witness += static_cast<int>(sizeof(controller) > 0);

    witness += ctrlpp::biquad<scalar>::low_pass(scalar{100}, scalar{1000}).has_value() ? 1 : 0;

    witness += ctrlpp::complementary_filter<scalar>::create(
                   ctrlpp::cf_config<scalar>{})
                   .has_value()
               ? 1
               : 0;

    witness += ctrlpp::online_planner_3rd<scalar>::create(
                   {.v_max = scalar{1}, .a_max = scalar{1}, .j_max = scalar{1}})
                   .has_value()
               ? 1
               : 0;

    return witness;
}
}

int main()
{
    return fold_all() > 0 ? 0 : 1;
}
PROBE

# --- Leg 1: host, exceptions and RTTI off -------------------------------------

echo "=== Leg 1: host no-exceptions (-fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS) ==="

"${host_cxx}" -std=c++20 -fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS \
    -I "${ctrlpp_include}" -isystem "${eigen_include}" \
    -c "${witness_tu}" -o /dev/null
echo "  float translation unit: OK"

"${host_cxx}" -std=c++20 -fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS \
    -I "${ctrlpp_include}" -isystem "${eigen_include}" \
    -c "${double_probe}" -o /dev/null
echo "  double probe:           OK"
echo "Leg 1 PASS"
echo

# --- Leg 2: host C++20, exceptions on (expected fallback) ---------------------

echo "=== Leg 2: host -std=c++20 exceptions on (ctrlpp::detail::expected fallback) ==="

"${host_cxx}" -std=c++20 \
    -I "${ctrlpp_include}" -isystem "${eigen_include}" \
    -c "${witness_tu}" -o /dev/null
echo "  float translation unit: OK"

"${host_cxx}" -std=c++20 \
    -I "${ctrlpp_include}" -isystem "${eigen_include}" \
    -c "${double_probe}" -o /dev/null
echo "  double probe:           OK"
echo "Leg 2 PASS"
echo

# --- Leg 3: arm-none-eabi Cortex-M7 cross-compile -----------------------------

echo "=== Leg 3: arm-none-eabi-g++ Cortex-M7 no-exceptions ==="

if ! command -v "${arm_cxx}" >/dev/null 2>&1; then
    echo "ERROR: ${arm_cxx} not found on PATH." >&2
    echo "  Install the arm-none-eabi GCC toolchain (expected at /usr/bin)." >&2
    exit 1
fi

# The bare-metal toolchain resolves <cmath> through its C library headers
# (newlib). Verify they are present so a missing runtime package produces a
# clear diagnostic instead of a cryptic "math.h: No such file" from cmath.
# A caller may point ARM_EXTRA_INCLUDE at a newlib include tree (added with
# -idirafter) when the headers live outside the toolchain sysroot.
arm_sysroot="$("${arm_cxx}" -print-sysroot)"
if [ -z "${ARM_EXTRA_INCLUDE:-}" ] && [ ! -f "${arm_sysroot}/include/math.h" ]; then
    echo "ERROR: the ${arm_cxx} C library headers (newlib) are not installed." >&2
    echo "  Expected ${arm_sysroot}/include/math.h." >&2
    echo "  Install the newlib package for this toolchain (arm-none-eabi-newlib)," >&2
    echo "  or set ARM_EXTRA_INCLUDE to a newlib include directory." >&2
    exit 1
fi

arm_extra_args=()
if [ -n "${ARM_EXTRA_INCLUDE:-}" ]; then
    arm_extra_args+=(-idirafter "${ARM_EXTRA_INCLUDE}")
    echo "  Using extra include: ${ARM_EXTRA_INCLUDE}"
fi

arm-none-eabi-g++ -std=c++20 -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS \
    ${arm_extra_args[@]+"${arm_extra_args[@]}"} \
    -I "${ctrlpp_include}" -isystem "${eigen_include}" \
    -c "${witness_tu}" -o /dev/null
echo "  float translation unit: OK"
echo "Leg 3 PASS"
echo

echo "All legs PASS (Leg 1 host no-exceptions, Leg 2 host C++20 fallback, Leg 3 arm-none-eabi Cortex-M7)."
