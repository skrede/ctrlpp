cmake_minimum_required(VERSION 3.28)
project(ctrlpp VERSION 0.3.6)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Applied to ctrlpp's own targets only (tests, compile tests, examples); the
# vendored test framework and the third-party backends are compiled with their
# own options and never receive these.
#
# The discarded-result diagnostic is promoted to an ERROR rather than left a
# warning: ctrlpp::expected is annotated at class level, so dropping it on the
# floor turns a reported failure into a silent no-op, and in a control loop that
# is a step which did nothing while the caller believes it succeeded. A warning
# is scrollback; an error is a gate. Both trees were measured to emit zero such
# diagnostics before the promotion, so nothing pre-existing is caught by it.
#
# /bigobj raises the object file's section-count ceiling from 65,279 to 2^32. Two
# Riccati test translation units crossed the lower ceiling on MSVC and stopped
# compiling at all. The switch sits on this shared block rather than on those two
# targets because the rest of the set is close behind, and that closeness was
# measured rather than assumed: compiling the same sources with g++ and
# per-function/per-data sections, and counting the COMDAT groups against the same
# 65,279, the two that failed emit 59,881 and 57,586, the heaviest one that still
# compiled emits 47,398, and three more sit between 41,799 and 47,398. A single
# added state dimension is worth roughly 7,700 groups, so the next unit to cross
# is one instantiation away and a per-target switch would have to be re-argued
# every time it happened. The g++ group count is a stand-in for the format's own
# count, not that count.
#
# What /bigobj does NOT fix is the compiler's peak resident set. Those same two
# units cost 4.55 GiB and 4.37 GiB to compile with no debug info, and this switch
# moves neither figure: it widens a field in the object file and does nothing
# else. The instantiation cost that made them heavy is untouched. Splitting the
# heavier one along its method tag was measured as an alternative and does not
# clear the ceiling either, because the piece that has to keep the shared cases
# still emits 55,125 groups -- an 8% reduction -- while the three pieces together
# cost 38% more compile time than the single unit does.
set(CTRLPP_WARNING_FLAGS
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive- /utf-8 /EHsc /we4834 /bigobj>
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
        -Werror=unused-result
        # gcc reports a discarded value of an attributed class TYPE under the promotion above, while clang reports it under the promotion below and reserves the one above for a discarded attributed FUNCTION result, so neither promotion alone enforces the gate on every supported compiler.
        -Werror=unused-value
    >
)

option(CTRLPP_BUILD_OSQP "Build OSQP solver backend for MPC" OFF)
option(CTRLPP_BUILD_NLOPT "Build NLopt solver backend for NMPC" OFF)
option(CTRLPP_BUILD_ARGMIN "Build with Argmin NLP solver backend" OFF)
set(CTRLPP_ARGMIN_SOURCE_DIR "" CACHE PATH "Path to local Argmin source checkout")
set(CTRLPP_ARGMIN_GIT_TAG "30a2cd4" CACHE STRING "Argmin git tag/branch/commit to fetch (default pinned to a literal SHA for reproducibility; this SHA selects the argmin revision that ships argmin/qp/, so the argmin_qp_solver linear-MPC QP backend is ACTIVE at this pin -- no local checkout needed. The library has no public users yet, so the default tracks current argmin development rather than lagging it; the bump path from the earlier fea57c2 pin is a clean ancestor-verified forward move, and this SHA carries the round-off-scaled QP infeasibility certificates plus the resolve hot-loop and polish-factor-reuse optimizations. It retains the compile-time constraint_count channel that keeps the bounded static NMPC path's result-multiplier storage inline. A named ref/tag replaces this development pin once an upstream release includes the required surface. Override with a branch name to track a tip)")
option(CTRLPP_CMAKE_FETCH_DEPS "Use CMake FetchContent to download and install dependencies" OFF)

option(CTRLPP_BUILD_TESTS "Build unit and compile tests" OFF)
option(CTRLPP_BUILD_FUZZ_TESTS "Build fuzz test targets" OFF)
option(CTRLPP_BUILD_PROPERTY_TESTS "Build property-based test targets (requires RapidCheck)" OFF)

# Selects which posture ctrlpp's own root test targets compile in. This is
# self-imposed dogfood of the throw-free discipline, not a consumer mandate: the
# library/interface/install targets never carry -fno-exceptions. OFF (default)
# builds ctrlpp's tests and Catch2 under -fno-exceptions -fno-rtti with
# CATCH_CONFIG_DISABLE_EXCEPTIONS. ON builds the exceptions carve-out: tests
# asserting that throwing convenience wrappers throw (they need REQUIRE_THROWS*,
# unavailable under disabled-exceptions Catch2) plus every OSQP/NLopt-linked test
# with Catch2 exceptions on. Benchmarks are a separate project whose embedded
# nanobench implementation currently requires exceptions.
option(CTRLPP_TESTS_WITH_EXCEPTIONS "Build ctrlpp's root tests in the exceptions carve-out instead of the default -fno-exceptions tree" OFF)

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
