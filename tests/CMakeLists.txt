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
        # An imported Catch2 binary cannot be recompiled with
        # CATCH_CONFIG_DISABLE_EXCEPTIONS, so the default no-exceptions tree cannot
        # honor its posture on it: -fno-exceptions test TUs would link an
        # exceptions-on Catch2 whose failure path throws and unwinds through
        # -fno-exceptions frames, terminating without a report. Require building
        # Catch2 from source instead.
        if (NOT CTRLPP_TESTS_WITH_EXCEPTIONS)
            message(FATAL_ERROR
                "The default no-exceptions test tree (CTRLPP_TESTS_WITH_EXCEPTIONS=OFF) "
                "needs Catch2 built from source with CATCH_CONFIG_DISABLE_EXCEPTIONS, but "
                "Catch2 was resolved via find_package. Reconfigure with "
                "-DCTRLPP_CMAKE_FETCH_DEPS=ON, or build the exceptions carve-out tree "
                "with -DCTRLPP_TESTS_WITH_EXCEPTIONS=ON.")
        endif ()
    endif ()

    # Match Catch2's exception posture to the tree so ctrlpp's own tests and their
    # test framework share one setting. The default tree compiles the fetched
    # Catch2 -fno-exceptions -fno-rtti with CATCH_CONFIG_DISABLE_EXCEPTIONS, so its
    # failure path reports-and-aborts (REQUIRE_THROWS* / CHECK_THROWS* become
    # unavailable); the exceptions carve-out tree leaves Catch2 exceptions-on. The
    # PUBLIC definition reaches every test TU through Catch2WithMain.
    if (CTRLPP_CMAKE_FETCH_DEPS AND NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
        target_compile_definitions(Catch2 PUBLIC CATCH_CONFIG_DISABLE_EXCEPTIONS)
        target_compile_options(Catch2 PRIVATE -fno-exceptions -fno-rtti)
        target_compile_options(Catch2WithMain PRIVATE -fno-exceptions -fno-rtti)
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
