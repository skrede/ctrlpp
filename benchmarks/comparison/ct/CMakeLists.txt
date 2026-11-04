# ct_optcon depends on ct_core and may need the kindr headers.
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

# Only GCC 13 compiles bench_nmpc_vs_ct: ct/optcon/dms/dms_core/TimeGrid.h:53 assigns a read-only
# member in a template body, which GCC 15 and 16 reject (-Wtemplate-body). Build it from a second
# configure with -DCMAKE_CXX_COMPILER=g++-13 and disclose that compiler with the row.
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
