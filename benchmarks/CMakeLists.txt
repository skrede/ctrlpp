cmake_minimum_required(VERSION 3.28)
project(ctrlpp_benchmarks)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 CONFIG REQUIRED)

# ctrlpp is header-only; create a lightweight interface target
add_library(ctrlpp_benchmarks_lib INTERFACE)
add_library(ctrlpp::ctrlpp ALIAS ctrlpp_benchmarks_lib)
target_include_directories(ctrlpp_benchmarks_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib/ctrlpp/include
)
target_link_libraries(ctrlpp_benchmarks_lib INTERFACE Eigen3::Eigen)

# nanobench microbenchmark framework
include(FetchContent)
FetchContent_Declare(nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG        v4.3.11
    GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(nanobench)

# Optional solver backends for MPC/MHE benchmarks
find_package(osqp QUIET)
find_package(nlopt QUIET)

if(osqp_FOUND)
    add_library(ctrlpp_benchmarks_osqp INTERFACE)
    add_library(ctrlpp::osqp ALIAS ctrlpp_benchmarks_osqp)
    target_link_libraries(ctrlpp_benchmarks_osqp INTERFACE ctrlpp::ctrlpp osqp::osqp)
    target_compile_definitions(ctrlpp_benchmarks_osqp INTERFACE CTRLPP_HAS_OSQP=1)
endif()

if(nlopt_FOUND)
    add_library(ctrlpp_benchmarks_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_benchmarks_nlopt)
    target_link_libraries(ctrlpp_benchmarks_nlopt INTERFACE ctrlpp::ctrlpp nlopt)
    target_compile_definitions(ctrlpp_benchmarks_nlopt INTERFACE CTRLPP_HAS_NLOPT=1)
endif()

add_subdirectory(internal)
add_subdirectory(comparison)
