include(${CMAKE_CURRENT_LIST_DIR}/../deps/argmin.cmake)

find_package(NLopt REQUIRED)

if(NOT TARGET ctrlpp::nlopt)
    add_library(ctrlpp_argmin_bench_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_argmin_bench_nlopt)
    target_link_libraries(ctrlpp_argmin_bench_nlopt INTERFACE ctrlpp::ctrlpp NLopt::nlopt)
    target_compile_definitions(ctrlpp_argmin_bench_nlopt INTERFACE CTRLPP_HAS_NLOPT=1)
endif()

set(ARGMIN_BENCH_TARGETS
    bench_slsqp
    bench_nmpc
    bench_nmhe
    bench_step_budget
    bench_profile
    bench_jacobian
    bench_schedule
    bench_sqp_variants
)

foreach(target IN LISTS ARGMIN_BENCH_TARGETS)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp")
        add_executable(${target} ${target}.cpp)
        target_link_libraries(${target} PRIVATE
            ctrlpp::ctrlpp
            ctrlpp::bench_problems
            ctrlpp::nlopt
            nanobench
            argmin::argmin
        )
    endif()
endforeach()

if(TARGET ctrlpp::osqp AND EXISTS "${argmin_SOURCE_DIR}/lib/argmin/include/argmin/qp/sparse_admm_qp.h")
    foreach(qp_bench bench_qp_vs_osqp bench_qp_vs_osqp_closedloop bench_qp_vs_osqp_churning bench_qp_periter)
        add_executable(${qp_bench} ${qp_bench}.cpp)
        target_link_libraries(${qp_bench} PRIVATE
            ctrlpp::ctrlpp
            ctrlpp::bench_problems
            ctrlpp::osqp
            nanobench
            argmin::argmin
        )
    endforeach()
else()
    message(STATUS "ctrlpp bench: skipping bench_qp_vs_osqp[_closedloop] (needs both OSQP and an argmin checkout shipping argmin/qp/sparse_admm_qp.h)")
endif()
