# qpOASES comparison benchmark.
#
# qpOASES is the canonical online active-set QP solver for embedded MPC.
# It is the de facto warm-started small-dense baseline.

if(NOT CTRLPP_FETCH_BENCHMARK_DEPS)
    message(FATAL_ERROR "The qpOASES comparison benchmark is fetched from upstream and must be opted into; re-run with -DCTRLPP_FETCH_BENCHMARK_DEPS=ON.")
endif()

include(FetchContent)

FetchContent_Declare(qpoases
    GIT_REPOSITORY https://github.com/coin-or/qpOASES.git
    GIT_TAG releases/3.2.2
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(QPOASES_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# qpOASES ships cmake_minimum_required(VERSION 2.6); lift to 3.5 so the
# configure step does not error on modern CMake.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)

FetchContent_MakeAvailable(qpoases)

add_executable(bench_qp_vs_qpoases bench_qp_vs_qpoases.cpp)

target_include_directories(bench_qp_vs_qpoases PRIVATE
    ${qpoases_SOURCE_DIR}/include
)

target_link_libraries(bench_qp_vs_qpoases PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    qpOASES
)
