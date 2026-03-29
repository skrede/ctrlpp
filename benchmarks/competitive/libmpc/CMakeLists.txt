# libmpc++ competitive benchmarks
# Dependency: find_package first, FetchContent fallback
find_package(mpc++ CONFIG QUIET)
if(NOT mpc++_FOUND)
    message(STATUS "libmpc++ not found via find_package, fetching via FetchContent")
    include(FetchContent)

    # Prevent libmpc++ from fetching its own OSQP if our target already exists
    if(NOT TARGET osqpstatic AND NOT TARGET osqp::osqp)
        set(LIBMPC_FETCH_OSQP ON)
    endif()

    FetchContent_Declare(libmpc
        GIT_REPOSITORY https://github.com/nicolapiccinelli/libmpc.git
        GIT_TAG        1.0.0
        GIT_SHALLOW    TRUE
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(libmpc)
endif()

add_executable(bench_mpc_vs_libmpc bench_mpc_vs_libmpc.cpp)
target_link_libraries(bench_mpc_vs_libmpc PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::osqp
    nanobench
    mpc++)
