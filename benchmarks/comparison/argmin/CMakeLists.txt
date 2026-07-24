# A local checkout (CTRLPP_ARGMIN_SOURCE_DIR) always wins; otherwise fetch only
# when opted in, else expect an installed argmin via find_package.
if(CTRLPP_ARGMIN_SOURCE_DIR)
    include(FetchContent)
    FetchContent_Declare(argmin
        SOURCE_DIR ${CTRLPP_ARGMIN_SOURCE_DIR}
        EXCLUDE_FROM_ALL)
    set(ARGMIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(argmin)
elseif(CTRLPP_FETCH_BENCHMARK_DEPS)
    include(FetchContent)
    FetchContent_Declare(argmin
        GIT_REPOSITORY https://github.com/skrede/argmin.git
        GIT_TAG ${CTRLPP_ARGMIN_GIT_TAG}
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL)
    set(ARGMIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(argmin)
else()
    find_package(argmin CONFIG REQUIRED)
endif()

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

# QP-vs-OSQP comparison: argmin's sparse ADMM QP against ctrlpp's OSQP backend.
# Unlike the SQP/NMPC benches above this needs OSQP, not NLopt, and it needs an
# argmin checkout that ships the QP headers. The default argmin pin now ships
# them, so a stock argmin build satisfies this; it is still guarded and skipped
# when either OSQP or the QP header is missing so a stock benchmark configure is
# unaffected.
if(TARGET ctrlpp::osqp AND EXISTS "${argmin_SOURCE_DIR}/lib/argmin/include/argmin/qp/sparse_admm_qp.h")
    foreach(qp_bench bench_qp_vs_osqp bench_qp_vs_osqp_closedloop bench_qp_vs_osqp_churning bench_qp_periter)
        add_executable(${qp_bench} ${qp_bench}.cpp)
        target_link_libraries(${qp_bench} PRIVATE
            ctrlpp::ctrlpp
            ctrlpp::osqp
            nanobench
            argmin::argmin
        )
    endforeach()
else()
    message(STATUS "ctrlpp bench: skipping bench_qp_vs_osqp[_closedloop] (needs both OSQP and an argmin checkout shipping argmin/qp/sparse_admm_qp.h)")
endif()
