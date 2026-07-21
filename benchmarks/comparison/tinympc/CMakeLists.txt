# TinyMPC comparison benchmark.
#
# TinyMPC is an embedded ADMM solver tuned for very small LMPC. Different
# design point from interior-point/active-set; the bench shows the tradeoff
# between embedded simplicity and convergence speed on a non-trivial problem.

if(NOT CTRLPP_FETCH_BENCHMARK_DEPS)
    message(FATAL_ERROR "The TinyMPC comparison benchmark is fetched from upstream and must be opted into; re-run with -DCTRLPP_FETCH_BENCHMARK_DEPS=ON.")
endif()

include(FetchContent)

FetchContent_Declare(tinympc
    GIT_REPOSITORY https://github.com/TinyMPC/TinyMPC.git
    GIT_TAG main
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

set(TINYMPC_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(TINYMPC_BUILD_TESTS    OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(tinympc)

add_executable(bench_lmpc_vs_tinympc bench_lmpc_vs_tinympc.cpp)

target_include_directories(bench_lmpc_vs_tinympc PRIVATE
    ${tinympc_SOURCE_DIR}/src
)

target_link_libraries(bench_lmpc_vs_tinympc PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    tinympcstatic
)
