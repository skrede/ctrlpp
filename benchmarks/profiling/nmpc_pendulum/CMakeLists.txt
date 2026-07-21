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
    add_library(ctrlpp_profiling_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_profiling_nlopt)
    target_link_libraries(ctrlpp_profiling_nlopt INTERFACE ctrlpp::ctrlpp NLopt::nlopt)
    target_compile_definitions(ctrlpp_profiling_nlopt INTERFACE CTRLPP_HAS_NLOPT=1)
endif()

set(PROFILE_TARGETS
    perf_nlopt_slsqp
    perf_argmin_slsqp
)

foreach(target IN LISTS PROFILE_TARGETS)
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} PRIVATE
        ctrlpp::ctrlpp
        ctrlpp::bench_problems
        ctrlpp::nlopt
        argmin::argmin
    )
endforeach()
