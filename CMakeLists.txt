cmake_minimum_required(VERSION 3.28)
project(ctrlpp VERSION 0.3.6)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CTRLPP_WARNING_FLAGS
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive- /utf-8 /EHsc>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
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
        -Wimplicit-fallthrough
        -Wformat=2
    >
)

option(CTRLPP_BUILD_OSQP "Build OSQP solver backend for MPC" OFF)
option(CTRLPP_BUILD_NLOPT "Build NLopt solver backend for NMPC" OFF)
option(CTRLPP_BUILD_ARGMIN "Build with Argmin NLP solver backend" OFF)
set(CTRLPP_ARGMIN_SOURCE_DIR "" CACHE PATH "Path to local Argmin source checkout")
set(CTRLPP_ARGMIN_GIT_TAG "30a2cd4" CACHE STRING "Argmin git tag/branch/commit to fetch (default pinned to a literal SHA for reproducibility; this SHA tracks the current tip of argmin milestone/v0.3.5 on origin and ships argmin/qp/, so the argmin_qp_solver linear-MPC QP backend is ACTIVE at this pin -- no local checkout needed. The library has no public users yet, so the default tracks argmin's tip rather than lagging it; the bump path from the earlier fea57c2 pin is a clean ancestor-verified forward move on the same branch, and this SHA carries the round-off-scaled QP infeasibility certificates plus the resolve hot-loop and polish-factor-reuse optimizations. It retains the compile-time constraint_count channel that keeps the bounded static NMPC path's result-multiplier storage inline. A named ref/tag replaces this floating SHA once argmin's milestone closes. Override with a branch name to track a tip)")
option(CTRLPP_CMAKE_FETCH_DEPS "Use CMake FetchContent to download and install dependencies" OFF)

option(CTRLPP_BUILD_TESTS "Build unit and compile tests" OFF)
option(CTRLPP_BUILD_FUZZ_TESTS "Build fuzz test targets" OFF)
option(CTRLPP_BUILD_PROPERTY_TESTS "Build property-based test targets (requires RapidCheck)" OFF)

# Selects which build tree ctrlpp's OWN tests and benches compile as. This is a
# self-imposed dogfood of the throw-free discipline, not a consumer mandate: the
# library/interface/install targets never carry -fno-exceptions. OFF (default)
# builds ctrlpp's test + bench targets AND Catch2 under -fno-exceptions -fno-rtti
# with CATCH_CONFIG_DISABLE_EXCEPTIONS. ON builds the exceptions carve-out --
# tests asserting the throwing convenience wrappers throw (they need
# REQUIRE_THROWS*, unavailable under disabled-exceptions Catch2) plus every
# OSQP/NLopt-linked test and bench -- with Catch2 exceptions-on.
option(CTRLPP_TESTS_WITH_EXCEPTIONS "Build ctrlpp's own tests/benches in the exceptions carve-out tree instead of the default -fno-exceptions tree" OFF)

option(CTRLPP_BUILD_EXAMPLES "Build examples" OFF)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_subdirectory(lib)

if (CTRLPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (CTRLPP_BUILD_TESTS OR CTRLPP_BUILD_FUZZ_TESTS OR CTRLPP_BUILD_PROPERTY_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
