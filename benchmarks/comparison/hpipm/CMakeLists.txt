# HPIPM comparison benchmark.
#
# Wires the HPIPM C API directly (the upstream hpipm-cpp wrapper does not
# build against current HPIPM struct definitions) and FetchContents both
# HPIPM and BLASFEO from upstream because the Arch hpipm-git / blasfeo-git
# packages ship LTO-only static archives with no usable symbol table.

if(NOT CTRLPP_FETCH_BENCHMARK_DEPS)
    message(FATAL_ERROR "The HPIPM comparison benchmark has no usable system package (upstream ships LTO-only archives) and must be fetched; re-run with -DCTRLPP_FETCH_BENCHMARK_DEPS=ON.")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/../deps/hpipm.cmake)

add_executable(bench_qp_vs_hpipm bench_qp_vs_hpipm.cpp)

target_include_directories(bench_qp_vs_hpipm PRIVATE
    ${blasfeo_SOURCE_DIR}/include
    ${hpipm_SOURCE_DIR}/include
)

target_link_libraries(bench_qp_vs_hpipm PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    hpipm
    blasfeo
)
