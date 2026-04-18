# ETH Control Toolbox comparison benchmarks
# Warning: ct_optcon depends on ct_core and may require kindr headers.
find_package(ct_optcon REQUIRED)

add_executable(bench_lqr_vs_ct bench_lqr_vs_ct.cpp)
target_link_libraries(bench_lqr_vs_ct PRIVATE
    ctrlpp::ctrlpp
    nanobench
    ct_optcon
)

add_executable(bench_dare_vs_ct bench_dare_vs_ct.cpp)
target_link_libraries(bench_dare_vs_ct PRIVATE
    ctrlpp::ctrlpp
    nanobench
    ct_optcon
)

add_executable(bench_care_vs_ct bench_care_vs_ct.cpp)
target_link_libraries(bench_care_vs_ct PRIVATE
    ctrlpp::ctrlpp
    nanobench
    ct_optcon
)

add_executable(bench_lqr_continuous_vs_ct bench_lqr_continuous_vs_ct.cpp)
target_link_libraries(bench_lqr_continuous_vs_ct PRIVATE
    ctrlpp::ctrlpp
    nanobench
    ct_optcon
)
