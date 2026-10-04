# qpOASES comparison benchmark.
#
# qpOASES is the canonical online active-set QP solver for embedded MPC.
# It is the de facto warm-started small-dense baseline.

if(NOT CTRLPP_FETCH_BENCHMARK_DEPS)
    message(FATAL_ERROR "The qpOASES comparison benchmark is fetched from upstream and must be opted into; re-run with -DCTRLPP_FETCH_BENCHMARK_DEPS=ON.")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/../deps/qpoases.cmake)

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
