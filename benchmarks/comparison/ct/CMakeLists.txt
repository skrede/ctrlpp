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

# bench_nmpc_vs_ct compiles only under GCC 13. The optcon umbrella pulls
# ct/optcon/dms/dms_core/TimeGrid.h, whose line 53 assigns the read-only member
# numberOfShots_ inside a template body; GCC 15 and 16 reject that with
# "error: assignment of read-only member 'numberOfShots_' [-Wtemplate-body]",
# while GCC 13 accepts it. Build the target from a second configure directory
# with -DCMAKE_CXX_COMPILER=g++-13 and disclose that compiler with the row. The
# target also needs the NLopt backend for the ctrlpp side of the comparison, so
# it is registered only where NLopt is installed.
find_package(NLopt QUIET)
if(TARGET NLopt::nlopt)
    add_executable(bench_nmpc_vs_ct bench_nmpc_vs_ct.cpp)
    target_link_libraries(bench_nmpc_vs_ct PRIVATE
        ctrlpp::ctrlpp
        nanobench
        ct_optcon
        NLopt::nlopt
    )
endif()
