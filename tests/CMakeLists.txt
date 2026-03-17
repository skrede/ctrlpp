if (CTRLPP_BUILD_TESTS)
    if (CTRLPP_CMAKE_FETCH_DEPS)
        FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.7.1
            EXCLUDE_FROM_ALL
        )
        FetchContent_MakeAvailable(Catch2)
    else ()
        find_package(Catch2 CONFIG REQUIRED)
    endif ()
    add_subdirectory(unit)
    add_subdirectory(compile)
    add_subdirectory(integration)
endif ()

if (CTRLPP_BUILD_FUZZ_TESTS)
    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR "Fuzz targets require Clang with libFuzzer support")
    endif ()
    add_subdirectory(fuzz)
endif ()
