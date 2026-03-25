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
    if (CTRLPP_BUILD_PROPERTY_TESTS)
        FetchContent_Declare(
            rapidcheck
            GIT_REPOSITORY https://github.com/emil-e/rapidcheck.git
            GIT_TAG b96a4e626ef4c7348dcd16c500353c2f997a9f3f
            EXCLUDE_FROM_ALL
        )
        FetchContent_MakeAvailable(rapidcheck)
        add_subdirectory(${rapidcheck_SOURCE_DIR}/extras/catch ${rapidcheck_BINARY_DIR}/extras/catch)

        # Suppress warnings from RapidCheck headers (third-party code)
        get_target_property(_rc_inc rapidcheck INTERFACE_INCLUDE_DIRECTORIES)
        if (_rc_inc)
            set_target_properties(rapidcheck PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_rc_inc}")
        endif ()
        get_target_property(_rc_catch_inc rapidcheck_catch INTERFACE_INCLUDE_DIRECTORIES)
        if (_rc_catch_inc)
            set_target_properties(rapidcheck_catch PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_rc_catch_inc}")
        endif ()
    endif ()
    add_subdirectory(unit)
    add_subdirectory(compile)
    if (NOT CTRLPP_CMAKE_FETCH_DEPS)
        add_subdirectory(integration)
    endif ()
endif ()

if (CTRLPP_BUILD_FUZZ_TESTS)
    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR "Fuzz targets require Clang with libFuzzer support")
    endif ()
    add_subdirectory(fuzz)
endif ()
