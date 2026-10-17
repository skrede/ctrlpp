# A local checkout (CTRLPP_ARGMIN_SOURCE_DIR) always wins; otherwise fetch only
# when opted in, else expect an installed argmin via find_package.
include(FetchContent)

if(CTRLPP_ARGMIN_SOURCE_DIR)
    FetchContent_Declare(argmin
        SOURCE_DIR ${CTRLPP_ARGMIN_SOURCE_DIR}
        EXCLUDE_FROM_ALL)
elseif(CTRLPP_FETCH_BENCHMARK_DEPS)
    set(_ctrlpp_bench_argmin_tag "${CTRLPP_ARGMIN_GIT_TAG}")
    set(_ctrlpp_bench_argmin_origin "passed in")
    if(NOT _ctrlpp_bench_argmin_tag)
        # The bench project cannot see the root project's cache, so the revision is read from the file
        # that declares it rather than restated; the configure re-runs when that file changes.
        set(_ctrlpp_root_lists "${CMAKE_CURRENT_LIST_DIR}/../../../CMakeLists.txt")
        set(_ctrlpp_pin_pattern "^set\\(CTRLPP_ARGMIN_GIT_TAG \"([^\"]+)\"")
        file(STRINGS "${_ctrlpp_root_lists}" _ctrlpp_pin_line REGEX "${_ctrlpp_pin_pattern}" LIMIT_COUNT 1)
        if(NOT _ctrlpp_pin_line MATCHES "${_ctrlpp_pin_pattern}")
            message(FATAL_ERROR "ctrlpp bench: no match for '${_ctrlpp_pin_pattern}' in ${_ctrlpp_root_lists}")
        endif()
        set(_ctrlpp_bench_argmin_tag "${CMAKE_MATCH_1}")
        set(_ctrlpp_bench_argmin_origin "the root project's pin")
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_ctrlpp_root_lists}")
    endif()
    message(STATUS "ctrlpp bench: fetching argmin ${_ctrlpp_bench_argmin_tag} (${_ctrlpp_bench_argmin_origin})")
    # A full clone, because a shallow one cannot always check out a bare commit SHA.
    FetchContent_Declare(argmin
        GIT_REPOSITORY https://github.com/skrede/argmin.git
        GIT_TAG ${_ctrlpp_bench_argmin_tag}
        EXCLUDE_FROM_ALL)
else()
    find_package(argmin CONFIG REQUIRED)
    return()
endif()

set(ARGMIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ARGMIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ARGMIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(argmin)
