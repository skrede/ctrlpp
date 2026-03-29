# ruckig competitive benchmarks
# Dependency: find_package first, FetchContent fallback
find_package(ruckig QUIET)
if(NOT ruckig_FOUND)
    message(STATUS "ruckig not found via find_package, fetching via FetchContent")
    include(FetchContent)
    FetchContent_Declare(ruckig
        GIT_REPOSITORY https://github.com/pantor/ruckig.git
        GIT_TAG        v0.17.3
        GIT_SHALLOW    TRUE
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(ruckig)
endif()

# Plan 03 adds benchmark executables here
