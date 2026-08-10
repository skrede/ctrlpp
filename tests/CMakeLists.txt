if (CTRLPP_BUILD_TESTS)
    # Prove with the active compiler that discarding ctrlpp::expected is an
    # error. A text search for the options cannot establish that the compiler
    # recognizes them or that the class-level attribute still emits the
    # diagnostic they promote.
    #
    # The non-MSVC branch carries both promotions because the compilers disagree
    # about which group a discarded value of an attributed class type belongs
    # to: clang files it under the second group and a discarded attributed
    # function result under the first, while gcc files both under the first.
    # Either promotion alone therefore leaves the gate off on some supported
    # compiler, so the pair is the gate and neither element is decorative.
    if (MSVC)
        set(CTRLPP_DISCARD_WARNING_ERROR_FLAG "/we4834")
    else ()
        set(CTRLPP_DISCARD_WARNING_ERROR_FLAG "-Werror=unused-result" "-Werror=unused-value")
    endif ()
    try_compile(
        CTRLPP_DISCARDED_EXPECTED_BASELINE_COMPILED
        SOURCES ${CMAKE_CURRENT_LIST_DIR}/compile/discarded_expected_canary.cpp
        CMAKE_FLAGS
            "-DCMAKE_CXX_STANDARD=20"
            "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
            "-DCMAKE_CXX_EXTENSIONS=OFF"
            "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/lib/ctrlpp/include"
        OUTPUT_VARIABLE CTRLPP_DISCARDED_EXPECTED_BASELINE_OUTPUT
    )
    if (NOT CTRLPP_DISCARDED_EXPECTED_BASELINE_COMPILED)
        message(FATAL_ERROR
            "The discarded-expected compile canary is not a valid positive "
            "control without warning promotion.\n"
            "${CTRLPP_DISCARDED_EXPECTED_BASELINE_OUTPUT}")
    endif ()
    try_compile(
        CTRLPP_DISCARDED_EXPECTED_COMPILED
        SOURCES ${CMAKE_CURRENT_LIST_DIR}/compile/discarded_expected_canary.cpp
        CMAKE_FLAGS
            "-DCMAKE_CXX_STANDARD=20"
            "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
            "-DCMAKE_CXX_EXTENSIONS=OFF"
            "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/lib/ctrlpp/include"
        COMPILE_DEFINITIONS ${CTRLPP_DISCARD_WARNING_ERROR_FLAG}
        OUTPUT_VARIABLE CTRLPP_DISCARDED_EXPECTED_OUTPUT
    )
    if (CTRLPP_DISCARDED_EXPECTED_COMPILED)
        message(FATAL_ERROR
            "Discarding ctrlpp::expected compiled successfully under the "
            "discard promotions ${CTRLPP_DISCARD_WARNING_ERROR_FLAG}; the "
            "discarded-result diagnostic is not enforced.\n"
            "${CTRLPP_DISCARDED_EXPECTED_OUTPUT}")
    endif ()

    # Every runtime test of the recursive ARX estimator passes whether or not the
    # estimator accepts input- and output-count template arguments, because no
    # legal instantiation ever supplied them. A refusal to compile is therefore
    # the only observable evidence that the unimplemented multi-channel surface
    # is gone.
    get_target_property(CTRLPP_CANARY_EIGEN_TARGET Eigen3::Eigen ALIASED_TARGET)
    if (NOT CTRLPP_CANARY_EIGEN_TARGET)
        set(CTRLPP_CANARY_EIGEN_TARGET Eigen3::Eigen)
    endif ()
    get_target_property(CTRLPP_CANARY_EIGEN_INCLUDE
        ${CTRLPP_CANARY_EIGEN_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
    try_compile(
        CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_COMPILED
        SOURCES ${CMAKE_CURRENT_LIST_DIR}/compile/recursive_arx_arity_canary.cpp
        CMAKE_FLAGS
            "-DCMAKE_CXX_STANDARD=20"
            "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
            "-DCMAKE_CXX_EXTENSIONS=OFF"
            "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/lib/ctrlpp/include;${CTRLPP_CANARY_EIGEN_INCLUDE}"
        OUTPUT_VARIABLE CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_OUTPUT
    )
    if (NOT CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_COMPILED)
        message(FATAL_ERROR
            "The recursive ARX arity canary is not a valid positive control: its "
            "three-argument arm does not compile.\n"
            "${CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_OUTPUT}")
    endif ()
    try_compile(
        CTRLPP_RECURSIVE_ARX_ARITY_COMPILED
        SOURCES ${CMAKE_CURRENT_LIST_DIR}/compile/recursive_arx_arity_canary.cpp
        CMAKE_FLAGS
            "-DCMAKE_CXX_STANDARD=20"
            "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
            "-DCMAKE_CXX_EXTENSIONS=OFF"
            "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/lib/ctrlpp/include;${CTRLPP_CANARY_EIGEN_INCLUDE}"
        COMPILE_DEFINITIONS -DCTRLPP_CANARY_FIVE_PARAMETER_ARITY
        OUTPUT_VARIABLE CTRLPP_RECURSIVE_ARX_ARITY_OUTPUT
    )
    if (CTRLPP_RECURSIVE_ARX_ARITY_COMPILED)
        message(FATAL_ERROR
            "A five-argument recursive ARX instantiation compiled successfully; "
            "the estimator still accepts the wider arity.\n"
            "${CTRLPP_RECURSIVE_ARX_ARITY_OUTPUT}")
    endif ()

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
    # Catch2 -fno-exceptions with CATCH_CONFIG_DISABLE_EXCEPTIONS, so its failure
    # path reports-and-aborts (REQUIRE_THROWS* / CHECK_THROWS* become
    # unavailable); the exceptions carve-out tree leaves Catch2 exceptions-on. The
    # PUBLIC definition reaches every test TU through Catch2WithMain.
    #
    # Catch2 keeps RTTI even though ctrlpp's own test TUs are compiled -fno-rtti.
    # RapidCheck's type-erased Any and its ShowType printer call typeid in headers
    # the property TUs instantiate, so those TUs cannot drop RTTI; and an
    # RTTI-enabled TU that decomposes an assertion emits the typeinfo of a
    # Catch::UnaryExpr/BinaryExpr instantiation, which references the typeinfo of
    # their base Catch::ITransientExpression. A -fno-rtti Catch2 never emits that
    # base typeinfo, so the reference goes unresolved at link time. Dropping the
    # flag on the framework costs nothing the discipline is about: ctrlpp's own
    # test sources still compile without exceptions and without RTTI.
    if (CTRLPP_CMAKE_FETCH_DEPS AND NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
        target_compile_definitions(Catch2 PUBLIC CATCH_CONFIG_DISABLE_EXCEPTIONS)
        target_compile_options(Catch2 PRIVATE -fno-exceptions)
        target_compile_options(Catch2WithMain PRIVATE -fno-exceptions)
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
    # The source purity gate and the self-test that proves it can fail. Both are
    # registered here at the top level rather than under integration/, which is
    # added only when dependencies are resolved by package discovery instead of
    # fetched, and would therefore be skipped in the configuration this project
    # builds by default.
    #
    # Neither carries a target dependency, unlike every other registration in
    # this tree: the gate reads source and consumes no build product. The
    # working directory is set because the script resolves its scan roots
    # relative to the repository root.
    #
    # Registration lives under the build-tests guard, so the gate is absent from
    # a library-only consumer build. That matches the intent -- it runs on every
    # test run, locally and in every workflow -- and is stated rather than
    # engineered around.
    # The interpreter is named explicitly rather than relying on the shebang,
    # because CTest launches the command directly and Windows cannot execute a
    # shell script that way. Where no interpreter exists the two tests are
    # skipped with a visible status line instead of registered and then failing
    # for a reason that has nothing to do with source purity. The gate checks
    # source text, which is platform-independent, so one platform running it
    # ratchets the property for all of them.
    find_program(CTRLPP_PURITY_SHELL bash)
    if (CTRLPP_PURITY_SHELL)
        add_test(
            NAME source_purity
            COMMAND ${CTRLPP_PURITY_SHELL}
                    ${CMAKE_SOURCE_DIR}/scripts/source_purity_check.sh
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
        add_test(
            NAME source_purity_self_test
            COMMAND ${CTRLPP_PURITY_SHELL}
                    ${CMAKE_SOURCE_DIR}/scripts/source_purity_check.sh --self-test
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
        set_tests_properties(source_purity source_purity_self_test PROPERTIES
            LABELS "purity"
        )
    else ()
        message(STATUS
            "No bash found: the source purity gate and its self-test are not registered.")
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
