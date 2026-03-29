# libmpc++ competitive benchmarks
# Dependency: find_package first, FetchContent fallback
find_package(mpc++ CONFIG QUIET)
if(NOT mpc++_FOUND)
    message(STATUS "libmpc++ not found via find_package, fetching via FetchContent")
    include(FetchContent)
    FetchContent_Declare(libmpc
        GIT_REPOSITORY https://github.com/nicolapiccinelli/libmpc.git
        GIT_TAG        1.0.0
        GIT_SHALLOW    TRUE
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(libmpc)
endif()

# Plan 03 adds benchmark executables here
