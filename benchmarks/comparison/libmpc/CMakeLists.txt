find_package(mpc++ CONFIG QUIET)

if (NOT mpc++_FOUND)
    message(STATUS "libmpc++ not found via find_package, fetching via FetchContent")
    include(FetchContent)

    if (NOT TARGET osqpstatic AND NOT TARGET osqp::osqp)
        set(LIBMPC_FETCH_OSQP ON)
    endif ()

    FetchContent_Declare(libmpc
        GIT_REPOSITORY https://github.com/nicolapiccinelli/libmpc.git
        GIT_TAG 1.0.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL)

    # libmpc++ has cmake_minimum_required(VERSION 3.0) which modern CMake rejects
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
    FetchContent_MakeAvailable(libmpc)
    unset(CMAKE_POLICY_VERSION_MINIMUM CACHE)

    # libmpc++ forgets to add its own include/ to the INTERFACE
    target_include_directories(mpc++ SYSTEM INTERFACE
        $<BUILD_INTERFACE:${libmpc_SOURCE_DIR}/include>)
endif ()

add_executable(bench_mpc_vs_libmpc bench_mpc_vs_libmpc.cpp)
target_link_libraries(bench_mpc_vs_libmpc PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    mpc++
)
