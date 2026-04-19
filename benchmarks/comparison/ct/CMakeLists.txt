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

# Note: a planned bench_nmpc_vs_ct (NLOC iLQR vs ctrlpp::nmpc) is not wired
# here because the installed ct_optcon 3.0.2 NLOC path does not compile on
# gcc 15: dms/TimeGrid.h assigns a const member, NLOptConSolver-impl.hpp
# overrides don't match the base virtuals under stricter checks, and the
# 2-arg constructor fails template instantiation. Revisit if ct_optcon ships
# a gcc 15-compatible release or a patched local build is available.
