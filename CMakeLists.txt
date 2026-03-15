cmake_minimum_required(VERSION 3.25)
project(ctrlpp VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CPPCTRL_WARNING_FLAGS
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive- /utf-8>
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

option(CPPCTRL_BUILD_TESTS "Build unit and compile tests" OFF)
option(CPPCTRL_BUILD_EXAMPLES "Build examples" OFF)
option(CPPCTRL_BUILD_FUZZ_TESTS "Build fuzz test targets" OFF)
option(CPPCTRL_CMAKE_FETCH_DEPS "Use CMake FetchContent to download and install dependencies" OFF)
option(CPPCTRL_ENABLE_EIGEN_POLICY "Build ctrlpp::eigen backed by Eigen" OFF)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_subdirectory(lib)

if (CPPCTRL_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (CPPCTRL_BUILD_TESTS OR CPPCTRL_BUILD_FUZZ_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
