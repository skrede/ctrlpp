# proxsuite's install path collides with FetchContent's EXCLUDE_FROM_ALL, so only an installed
# proxsuite is used.

find_package(proxsuite REQUIRED)

add_executable(bench_qp_vs_proxqp bench_qp_vs_proxqp.cpp)

target_link_libraries(bench_qp_vs_proxqp PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    proxsuite::proxsuite
)
