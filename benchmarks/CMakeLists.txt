cmake_minimum_required(VERSION 3.28)
project(ctrlpp_benchmarks)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CTRLPP_FETCH_BENCHMARK_DEPS "Fetch missing benchmark dependencies via FetchContent instead of only using installed ones" OFF)

# The root never adds this project, so it mirrors the unit trees' switch. nanobench, embedded in
# every internal bench TU, has no disable-exceptions mode, so the benches build exceptions-on either way.
option(CTRLPP_TESTS_WITH_EXCEPTIONS "Mirror of the unit trees' switch for the standalone bench project (see note: nanobench currently forces exceptions-on)" OFF)

find_package(Eigen3 CONFIG REQUIRED)

add_library(ctrlpp_benchmarks_lib INTERFACE)
add_library(ctrlpp::ctrlpp ALIAS ctrlpp_benchmarks_lib)
# The benchmark root is on the path too, so every harness can reach the shared
# bench_construct.h without a relative include.
target_include_directories(ctrlpp_benchmarks_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib/ctrlpp/include
)
target_link_libraries(ctrlpp_benchmarks_lib INTERFACE Eigen3::Eigen)

add_library(ctrlpp_bench_problems INTERFACE)
add_library(ctrlpp::bench_problems ALIAS ctrlpp_bench_problems)
target_include_directories(ctrlpp_bench_problems INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/problems
)
target_link_libraries(ctrlpp_bench_problems INTERFACE ctrlpp::ctrlpp)

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

# The cross-solver QP suite builds OSQP from source so that every solver it compares shares one set
# of compiler flags; every other configuration keeps the installed package, which osqp-eigen also links.
if(CTRLPP_FETCH_BENCHMARK_DEPS AND CTRLPP_BENCH_QP_SOLVERS)
    include(${CMAKE_CURRENT_SOURCE_DIR}/comparison/deps/osqp.cmake)
    set(ctrlpp_bench_osqp_target osqp::osqpstatic)
else()
    find_package(osqp QUIET)
    if(osqp_FOUND)
        set(ctrlpp_bench_osqp_target osqp::osqp)
    endif()
endif()
find_package(NLopt QUIET)

if(ctrlpp_bench_osqp_target)
    add_library(ctrlpp_benchmarks_osqp INTERFACE)
    add_library(ctrlpp::osqp ALIAS ctrlpp_benchmarks_osqp)
    target_link_libraries(ctrlpp_benchmarks_osqp INTERFACE ctrlpp::ctrlpp ${ctrlpp_bench_osqp_target})
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
