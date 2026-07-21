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
