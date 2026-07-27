cmake_minimum_required(VERSION 3.28)
project(ctrlpp_benchmarks)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Benchmark dependencies are resolved with find_package by default (assume the
# dependency is installed on the system). Set CTRLPP_FETCH_BENCHMARK_DEPS=ON to
# instead download and build the missing ones via FetchContent.
option(CTRLPP_FETCH_BENCHMARK_DEPS "Fetch missing benchmark dependencies via FetchContent instead of only using installed ones" OFF)

# This benchmarks project is standalone (the root CMakeLists never adds it), so
# the root CTRLPP_TESTS_WITH_EXCEPTIONS is invisible here and the bench tree
# carries its own switch mirroring the unit trees. The ctrlpp bench logic is
# throw-free, but every internal bench TU embeds nanobench's implementation, and
# nanobench has no disable-exceptions mode, so -fno-exceptions cannot currently be
# honored on the benches (see benchmarks/internal/CMakeLists.txt). The switch is
# kept for parity and for when a no-exceptions microbench path exists; the
# comparison/profiling benches link throwing solvers and stay exceptions-on.
option(CTRLPP_TESTS_WITH_EXCEPTIONS "Mirror of the unit trees' switch for the standalone bench project (see note: nanobench currently forces exceptions-on)" OFF)

find_package(Eigen3 CONFIG REQUIRED)

# ctrlpp is header-only; create a lightweight interface target
add_library(ctrlpp_benchmarks_lib INTERFACE)
add_library(ctrlpp::ctrlpp ALIAS ctrlpp_benchmarks_lib)
# The benchmark root is on the path too, so every harness can reach the shared
# bench_construct.h without a relative include.
target_include_directories(ctrlpp_benchmarks_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib/ctrlpp/include
)
target_link_libraries(ctrlpp_benchmarks_lib INTERFACE Eigen3::Eigen)

# Shared problem corpus consumed by both comparison/ and profiling/ harnesses.
# Header-only; do not put cpp files here.
add_library(ctrlpp_bench_problems INTERFACE)
add_library(ctrlpp::bench_problems ALIAS ctrlpp_bench_problems)
target_include_directories(ctrlpp_bench_problems INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/problems
)
target_link_libraries(ctrlpp_bench_problems INTERFACE ctrlpp::ctrlpp)

# nanobench microbenchmark framework
if(CTRLPP_FETCH_BENCHMARK_DEPS)
    include(FetchContent)
    FetchContent_Declare(nanobench
        GIT_REPOSITORY https://github.com/martinus/nanobench.git
        GIT_TAG        v4.3.11
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(nanobench)
else()
    find_package(nanobench CONFIG REQUIRED)
endif()

# Optional solver backends for MPC/MHE benchmarks
find_package(osqp QUIET)
find_package(NLopt QUIET)

if(osqp_FOUND)
    add_library(ctrlpp_benchmarks_osqp INTERFACE)
    add_library(ctrlpp::osqp ALIAS ctrlpp_benchmarks_osqp)
    target_link_libraries(ctrlpp_benchmarks_osqp INTERFACE ctrlpp::ctrlpp osqp::osqp)
    target_compile_definitions(ctrlpp_benchmarks_osqp INTERFACE CTRLPP_HAS_OSQP=1)
endif()

if(NLopt_FOUND)
    add_library(ctrlpp_benchmarks_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_benchmarks_nlopt)
    target_link_libraries(ctrlpp_benchmarks_nlopt INTERFACE ctrlpp::ctrlpp NLopt::nlopt)
    target_compile_definitions(ctrlpp_benchmarks_nlopt INTERFACE CTRLPP_HAS_NLOPT=1)
endif()

add_subdirectory(internal)
add_subdirectory(comparison)
add_subdirectory(profiling)
