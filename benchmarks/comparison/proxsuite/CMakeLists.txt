# proxsuite (ProxQP) comparison benchmark.
#
# proxsuite is the modern proximal QP solver from Simple Robotics. The dense
# variant is appropriate for small MPC-shaped QPs and is the closest
# apples-to-apples comparison vs OSQP/ctrlpp's wrapper.
#
# proxsuite's CMake install path collides with FetchContent's EXCLUDE_FROM_ALL
# semantics, so prefer the system install. Arch package: `proxsuite`.

find_package(proxsuite REQUIRED)

add_executable(bench_qp_vs_proxqp bench_qp_vs_proxqp.cpp)

target_link_libraries(bench_qp_vs_proxqp PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    proxsuite::proxsuite
)
