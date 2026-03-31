find_package(OsqpEigen REQUIRED)

add_executable(bench_qp_vs_osqp_eigen bench_qp_vs_osqp_eigen.cpp)
target_link_libraries(bench_qp_vs_osqp_eigen PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::osqp
    nanobench
    OsqpEigen::OsqpEigen
)
