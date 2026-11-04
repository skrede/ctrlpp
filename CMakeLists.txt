cmake_minimum_required(VERSION 3.28)
project(ctrlpp VERSION 0.3.6)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Applied to ctrlpp's own targets only, never to the vendored test framework or the backends. The
# discard promotions make a discarded standard-library [[nodiscard]] result, and an expression
# statement whose value is unused, an error on ctrlpp's targets.
set(CTRLPP_WARNING_FLAGS
    # /bigobj lifts MSVC's 65,279-section object ceiling, which two Riccati test units crossed. It is
    # on the shared set because the other units sit about one instantiation below the ceiling, and
    # splitting a failing unit was measured not to clear it.
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
        -Werror=unused-value
    >
    # gcc reports provably dead one-row branches of Eigen 3.4's Householder solves
    # (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=106247; Eigen issues 2506 and 2906, fixed on Eigen master by f7772e3)
    # and LDLT solves (unfixed, GCC PR56456) under either flag. Tree-wide because the need follows from instantiated
    # dimensions the build cannot see; clang's -Wsometimes-uninitialized and -Warray-bounds remain the backstop.
    $<$<CXX_COMPILER_ID:GNU>:-Wno-array-bounds -Wno-maybe-uninitialized>
)

option(CTRLPP_BUILD_OSQP "Build OSQP solver backend for MPC" OFF)
option(CTRLPP_BUILD_NLOPT "Build NLopt solver backend for NMPC" OFF)
option(CTRLPP_BUILD_ARGMIN "Build with Argmin NLP solver backend" OFF)
set(CTRLPP_ARGMIN_SOURCE_DIR "" CACHE PATH "Path to local Argmin source checkout")
set(CTRLPP_ARGMIN_GIT_TAG "0a638cb87b6e959f7cb7b49f0afa34f137b93b69" CACHE STRING "Argmin git tag/branch/commit to fetch. The default is a full commit SHA for reproducibility. It selects argmin 0.3.5, a forward move from the scalar-first template revision whose ancestry was verified, which carries the argmin/qp/ sparse operator-splitting QP backend with a 32-bit, domain-validated iteration budget and a nothrow-movable solver, so linear MPC and MHE hold argmin_qp_solver by value, an active-set equality subproblem that solves without heap allocation below 48 variables (larger subproblems still allocate in the periodic full QR refactorization), a multiplier estimate and an equality restoration that no longer allocate (other per-iteration allocations in the SQP step remain), the allocation-free at-solution null step of the NW-SQP policies, and zero-initialized L-BFGS-B workspaces. A named release tag replaces this development pin once an upstream release carries this surface. Override with a branch name or commit SHA to track another revision")
# An installed argmin older than this predates the scalar-first templates and the movable sparse QP solver the adapters use.
set(CTRLPP_ARGMIN_MIN_VERSION "0.3.5")
option(CTRLPP_CMAKE_FETCH_DEPS "Use CMake FetchContent to download and install dependencies" OFF)

option(CTRLPP_BUILD_TESTS "Build unit and compile tests" OFF)
option(CTRLPP_BUILD_FUZZ_TESTS "Build fuzz test targets" OFF)
option(CTRLPP_BUILD_PROPERTY_TESTS "Build property-based test targets (requires RapidCheck and CTRLPP_TESTS_WITH_EXCEPTIONS)" OFF)

# The -fno-exceptions -fno-rtti default is self-imposed on ctrlpp's own tests, not a consumer
# mandate: the library and install targets never carry it. ON builds the tests that need exceptions
# (REQUIRE_THROWS*, unavailable in disabled-exceptions Catch2) and every OSQP/NLopt-linked test.
option(CTRLPP_TESTS_WITH_EXCEPTIONS "Build ctrlpp's root tests in the exceptions carve-out instead of the default -fno-exceptions tree" OFF)

# Scoped to the three single-precision targets because the volume of these diagnostics elsewhere has
# never been measured. clang needs a third flag: gcc's float-conversion also covers the narrowing
# into single precision, which clang keeps in a group of its own.
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

# Not default-OFF like the build options: a standalone build is usually installed, while a
# subproject's default keeps a parent's prefix to what the parent asked for. A parent overrides it
# by setting the variable before the add_subdirectory that adds ctrlpp; afterwards is too late.
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
