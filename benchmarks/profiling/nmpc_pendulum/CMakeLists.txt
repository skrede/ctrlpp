include(${CMAKE_CURRENT_LIST_DIR}/../../comparison/deps/argmin.cmake)

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
