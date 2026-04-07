include(FetchContent)

if(CTRLPP_ARGMIN_SOURCE_DIR)
    FetchContent_Declare(nablapp
        SOURCE_DIR ${CTRLPP_ARGMIN_SOURCE_DIR}
        EXCLUDE_FROM_ALL)
else()
    FetchContent_Declare(nablapp
        GIT_REPOSITORY https://github.com/skrede/nablapp.git
        GIT_TAG milestone/v0.2.1
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL)
endif()

set(NABLAPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(NABLAPP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(NABLAPP_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(nablapp)

find_package(NLopt REQUIRED)

if(NOT TARGET ctrlpp::nlopt)
    add_library(ctrlpp_argmin_bench_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_argmin_bench_nlopt)
    target_link_libraries(ctrlpp_argmin_bench_nlopt INTERFACE ctrlpp::ctrlpp nlopt)
    target_compile_definitions(ctrlpp_argmin_bench_nlopt INTERFACE CTRLPP_HAS_NLOPT=1)
endif()

set(ARGMIN_BENCH_TARGETS
    bench_slsqp
    bench_mma
    bench_cobyla
    bench_isres
    bench_nmpc
    bench_nmhe
    bench_step_budget
    bench_lbfgsb
    bench_profile
    bench_jacobian
    bench_schedule
)

foreach(target IN LISTS ARGMIN_BENCH_TARGETS)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp")
        add_executable(${target} ${target}.cpp)
        target_link_libraries(${target} PRIVATE
            ctrlpp::ctrlpp
            ctrlpp::nlopt
            nanobench
            nablapp::nablapp
        )
    endif()
endforeach()
