# ETH Control Toolbox competitive benchmarks
# No FetchContent fallback -- catkin heritage, AUR available (yay -S control-toolbox)
# Warning: ct_optcon depends on ct_core and may require kindr headers.
find_package(ct_optcon REQUIRED)

add_executable(bench_lqr_vs_ct bench_lqr_vs_ct.cpp)
target_link_libraries(bench_lqr_vs_ct PRIVATE
    ctrlpp::ctrlpp
    nanobench
    ct_optcon)
