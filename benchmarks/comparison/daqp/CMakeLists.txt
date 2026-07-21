# DAQP comparison benchmark.
#
# DAQP is a dual active-set QP solver designed for embedded MPC inner subproblems.

if(NOT CTRLPP_FETCH_BENCHMARK_DEPS)
    message(FATAL_ERROR "The DAQP comparison benchmark is fetched from upstream and must be opted into; re-run with -DCTRLPP_FETCH_BENCHMARK_DEPS=ON.")
endif()

include(FetchContent)

FetchContent_Declare(daqp
    GIT_REPOSITORY https://github.com/darnstrom/daqp.git
    GIT_TAG v0.7.2
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

set(DAQP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(DAQP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(daqp)

add_executable(bench_qp_vs_daqp bench_qp_vs_daqp.cpp)

target_include_directories(bench_qp_vs_daqp PRIVATE
    ${daqp_SOURCE_DIR}/include
)

target_link_libraries(bench_qp_vs_daqp PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    daqpstat
)
