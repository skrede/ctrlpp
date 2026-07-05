cmake_minimum_required(VERSION 3.28)
project(ctrlpp VERSION 0.3.6)

set(CMAKE_CXX_STANDARD 23)
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
set(CTRLPP_ARGMIN_GIT_TAG "6219d0059fdbe3699cacd7f8e9fae4e8b720fd5a" CACHE STRING "Argmin git tag/branch/commit to fetch (default pinned to literal SHA for reproducibility; override with milestone/v0.3.0 to track branch tip)")
option(CTRLPP_CMAKE_FETCH_DEPS "Use CMake FetchContent to download and install dependencies" OFF)

option(CTRLPP_BUILD_TESTS "Build unit and compile tests" OFF)
option(CTRLPP_BUILD_FUZZ_TESTS "Build fuzz test targets" OFF)
option(CTRLPP_BUILD_PROPERTY_TESTS "Build property-based test targets (requires RapidCheck)" OFF)

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
