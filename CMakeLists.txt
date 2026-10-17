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
    # gcc 15 and 16 warn on provably dead branches of Eigen 3.4's Householder and LDLT solves at one-row operands
    # (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=106247 under PR56456; Eigen issues 2506 and 2906, fixed on Eigen master
    # by f7772e3); silencing one of the two reports the same code under the other. Remove once the fetched Eigen is clean.
    $<$<CXX_COMPILER_ID:GNU>:-Wno-array-bounds -Wno-maybe-uninitialized>
)

option(CTRLPP_BUILD_OSQP "Build OSQP solver backend for MPC" OFF)
option(CTRLPP_BUILD_NLOPT "Build NLopt solver backend for NMPC" OFF)
option(CTRLPP_BUILD_ARGMIN "Build with Argmin NLP solver backend" OFF)
set(CTRLPP_ARGMIN_SOURCE_DIR "" CACHE PATH "Path to local Argmin source checkout")
set(CTRLPP_ARGMIN_GIT_TAG "fcfb88b477899ed905c9504293bccc0020de1cc6" CACHE STRING "Argmin git tag/branch/commit to fetch. The default is a full commit SHA for reproducibility. It selects argmin's scalar-first template revision, a forward move whose ancestry from the previous default was verified, which carries the argmin/qp/ sparse operator-splitting QP backend with a 32-bit, domain-validated iteration budget, the allocation-free at-solution null step of the NW-SQP policies, and zero-initialized L-BFGS-B workspaces. A named release tag replaces this development pin once an upstream release carries this surface. Override with a branch name or commit SHA to track another revision")
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

# Raises the promotion and narrowing diagnostics from warnings to errors on the
# three targets that instantiate the library at single precision, and on nothing
# else. The promotion is target-scoped rather than global because the volume
# outside that tier has never been measured: promoting across the whole tree
# would turn an unmeasured population of pre-existing diagnostics into build
# failures. Default OFF, so an ordinary build sees exactly what it saw before.
#
# The clang list carries a third promotion because the compilers group these
# diagnostics differently: gcc's float-conversion covers both the narrowing into
# single precision and a float-to-integer conversion, while clang keeps the
# narrowing in a group gcc has no option for at all.
option(CTRLPP_FLOAT_TIER_STRICT "Promote the promotion and narrowing diagnostics to errors on the single-precision targets" OFF)
set(CTRLPP_FLOAT_TIER_STRICT_FLAGS "")
if (CTRLPP_FLOAT_TIER_STRICT AND NOT MSVC)
    set(CTRLPP_FLOAT_TIER_STRICT_FLAGS -Werror=double-promotion -Werror=float-conversion)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        list(APPEND CTRLPP_FLOAT_TIER_STRICT_FLAGS -Werror=implicit-float-conversion)
    endif ()
endif ()

option(CTRLPP_BUILD_EXAMPLES "Build examples" OFF)

option(CTRLPP_BUILD_VALIDATION "Build the Octave cross-validation harness" OFF)

# Deliberately not following the default-OFF convention the CTRLPP_BUILD_* options above
# share: those decide what this build contains, and this decides where it ends up. Standalone
# builds are asked to be installed often enough that OFF would be a trap; a parent that adds
# ctrlpp as a subproject gets OFF, so a consumer's prefix carries only what the consumer
# asked for. A parent overrides the default by setting the variable before the
# add_subdirectory that adds ctrlpp; setting it afterwards is too late.
option(CTRLPP_INSTALL "Install ctrlpp into the install prefix" ${PROJECT_IS_TOP_LEVEL})

# Config-package gate flags. lib/ctrlpp/ raises each to TRUE when the backend joins the
# export set, and the config file template reads them at substitution time.
set(CTRLPP_OSQP_INSTALLED FALSE)
set(CTRLPP_NLOPT_INSTALLED FALSE)
set(CTRLPP_ARGMIN_INSTALLED FALSE)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_subdirectory(lib)

if (CTRLPP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (CTRLPP_BUILD_TESTS OR CTRLPP_BUILD_FUZZ_TESTS OR CTRLPP_BUILD_PROPERTY_TESTS OR CTRLPP_BUILD_VALIDATION)
    enable_testing()
endif ()

if (CTRLPP_BUILD_TESTS OR CTRLPP_BUILD_FUZZ_TESTS OR CTRLPP_BUILD_PROPERTY_TESTS)
    add_subdirectory(tests)
endif ()

# After lib, because the cross-validation cases link the library target this
# project builds rather than one the harness declares for itself.
if (CTRLPP_BUILD_VALIDATION)
    add_subdirectory(validation)
endif ()
