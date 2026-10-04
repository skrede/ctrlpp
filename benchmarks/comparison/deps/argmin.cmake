# A local checkout (CTRLPP_ARGMIN_SOURCE_DIR) always wins; otherwise fetch only
# when opted in, else expect an installed argmin via find_package.
include(FetchContent)

if(CTRLPP_ARGMIN_SOURCE_DIR)
    FetchContent_Declare(argmin
        SOURCE_DIR ${CTRLPP_ARGMIN_SOURCE_DIR}
        EXCLUDE_FROM_ALL)
elseif(CTRLPP_FETCH_BENCHMARK_DEPS)
    # A full clone, because a shallow one cannot always check out a bare commit SHA.
    FetchContent_Declare(argmin
        GIT_REPOSITORY https://github.com/skrede/argmin.git
        GIT_TAG ${CTRLPP_ARGMIN_GIT_TAG}
        EXCLUDE_FROM_ALL)
else()
    find_package(argmin CONFIG REQUIRED)
    return()
endif()

set(ARGMIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ARGMIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ARGMIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(argmin)
