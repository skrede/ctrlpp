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

    # Prove with the active compiler that a single-to-double promotion and a
    # narrowing into single precision are errors under the promotions the
    # single-precision targets carry when CTRLPP_FLOAT_TIER_STRICT is on. A text
    # search for the options cannot establish that the compiler recognizes them
    # or that it still emits the diagnostics they promote, and the tier is clean,
    # so a green build says nothing on its own.
    #
    # The clang branch carries a third promotion because the two compilers group
    # these diagnostics differently: gcc's float-conversion covers both the
    # narrowing into single precision and a float-to-integer conversion, while
    # clang keeps the narrowing in a group gcc has no option for at all -- so one
    # shared list either misses the narrowing on clang or fails to parse on gcc.
    # The Microsoft compiler promotes neither diagnostic on its own, and a pair
    # run with no promotion to apply would report the absence of an option as an
    # unenforced gate.
    if (MSVC)
        message(STATUS
            "ctrlpp: the Microsoft compiler has no single-warning promotion for the "
            "single-precision diagnostics, skipping the float promotion canary pair.")
    else ()
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set(CTRLPP_FLOAT_WARNING_ERROR_FLAG "-Werror=double-promotion"
                "-Werror=float-conversion" "-Werror=implicit-float-conversion")
        else ()
            set(CTRLPP_FLOAT_WARNING_ERROR_FLAG "-Werror=double-promotion"
                "-Werror=float-conversion")
        endif ()
        try_compile(
            CTRLPP_FLOAT_PROMOTION_BASELINE_COMPILED
            SOURCES ${CMAKE_CURRENT_LIST_DIR}/compile/float_promotion_canary.cpp
            CMAKE_FLAGS
                "-DCMAKE_CXX_STANDARD=20"
                "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
                "-DCMAKE_CXX_EXTENSIONS=OFF"
            OUTPUT_VARIABLE CTRLPP_FLOAT_PROMOTION_BASELINE_OUTPUT
        )
        if (NOT CTRLPP_FLOAT_PROMOTION_BASELINE_COMPILED)
            message(FATAL_ERROR
                "The float promotion compile canary is not a valid positive "
                "control without warning promotion.\n"
                "${CTRLPP_FLOAT_PROMOTION_BASELINE_OUTPUT}")
        endif ()
        try_compile(
            CTRLPP_FLOAT_PROMOTION_COMPILED
            SOURCES ${CMAKE_CURRENT_LIST_DIR}/compile/float_promotion_canary.cpp
            CMAKE_FLAGS
                "-DCMAKE_CXX_STANDARD=20"
                "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
                "-DCMAKE_CXX_EXTENSIONS=OFF"
            COMPILE_DEFINITIONS ${CTRLPP_FLOAT_WARNING_ERROR_FLAG}
            OUTPUT_VARIABLE CTRLPP_FLOAT_PROMOTION_OUTPUT
        )
        if (CTRLPP_FLOAT_PROMOTION_COMPILED)
            message(FATAL_ERROR
                "A single-to-double promotion and a narrowing into single "
                "precision compiled successfully under the promotions "
                "${CTRLPP_FLOAT_WARNING_ERROR_FLAG}; the single-precision "
                "diagnostics are not enforced.\n"
                "${CTRLPP_FLOAT_PROMOTION_OUTPUT}")
        endif ()
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
    function(ctrlpp_try_canary result source)
        try_compile(${result}
            SOURCES ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/compile/${source}
            CMAKE_FLAGS
                "-DCMAKE_CXX_STANDARD=20"
                "-DCMAKE_CXX_STANDARD_REQUIRED=ON"
                "-DCMAKE_CXX_EXTENSIONS=OFF"
                "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/lib/ctrlpp/include;${CTRLPP_CANARY_EIGEN_INCLUDE}"
            COMPILE_DEFINITIONS ${ARGN}
            OUTPUT_VARIABLE output)
        set(${result} ${${result}} PARENT_SCOPE)
        set(${result}_OUTPUT "${output}" PARENT_SCOPE)
    endfunction()
    ctrlpp_try_canary(CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_COMPILED recursive_arx_arity_canary.cpp)
    if (NOT CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_COMPILED)
        message(FATAL_ERROR
            "The recursive ARX arity canary is not a valid positive control: its "
            "three-argument arm does not compile.\n"
            "${CTRLPP_RECURSIVE_ARX_ARITY_BASELINE_COMPILED_OUTPUT}")
    endif ()
    ctrlpp_try_canary(CTRLPP_RECURSIVE_ARX_ARITY_COMPILED recursive_arx_arity_canary.cpp -DCTRLPP_CANARY_FIVE_PARAMETER_ARITY)
    if (CTRLPP_RECURSIVE_ARX_ARITY_COMPILED)
        message(FATAL_ERROR
            "A five-argument recursive ARX instantiation compiled successfully; "
            "the estimator still accepts the wider arity.\n"
            "${CTRLPP_RECURSIVE_ARX_ARITY_COMPILED_OUTPUT}")
    endif ()

    # No runtime test can observe a constructor that is never called, so a
    # refusal to compile is the only observable evidence that the unchecked PID
    # construction path is closed.
    ctrlpp_try_canary(CTRLPP_PID_CONSTRUCTOR_BASELINE_COMPILED pid_private_constructor_canary.cpp)
    if (NOT CTRLPP_PID_CONSTRUCTOR_BASELINE_COMPILED)
        message(FATAL_ERROR
            "The pid constructor canary is not a valid positive control: its "
            "create arm does not compile.\n"
            "${CTRLPP_PID_CONSTRUCTOR_BASELINE_COMPILED_OUTPUT}")
    endif ()
    ctrlpp_try_canary(CTRLPP_PID_DIRECT_CONSTRUCTION_COMPILED pid_private_constructor_canary.cpp -DCTRLPP_CANARY_PID_DIRECT_CONSTRUCTION)
    if (CTRLPP_PID_DIRECT_CONSTRUCTION_COMPILED)
        message(FATAL_ERROR
            "A direct pid construction compiled; the constructor is not private.\n"
            "${CTRLPP_PID_DIRECT_CONSTRUCTION_COMPILED_OUTPUT}")
    endif ()

    # The particle filter's constructor is public so the factory can build the filter in the caller's storage,
    # which leaves the passkey as its only guard: a direct construction must fail, and fail on the tag.
    ctrlpp_try_canary(CTRLPP_PF_PASSKEY_BASELINE_COMPILED particle_filter_passkey_canary.cpp)
    if (NOT CTRLPP_PF_PASSKEY_BASELINE_COMPILED)
        message(FATAL_ERROR
            "The particle filter passkey canary is not a valid positive control: "
            "its create arm does not compile.\n"
            "${CTRLPP_PF_PASSKEY_BASELINE_COMPILED_OUTPUT}")
    endif ()
    ctrlpp_try_canary(CTRLPP_PF_DIRECT_CONSTRUCTION_COMPILED particle_filter_passkey_canary.cpp -DCTRLPP_CANARY_PF_DIRECT_CONSTRUCTION)
    set(output "${CTRLPP_PF_DIRECT_CONSTRUCTION_COMPILED_OUTPUT}")
    if (CTRLPP_PF_DIRECT_CONSTRUCTION_COMPILED)
        message(FATAL_ERROR "A direct particle filter construction compiled; the passkey is not enforced.\n${output}")
    elseif (NOT output MATCHES "validated_tag::validated_tag|constructor of class '[^']*::validated_tag'")
        message(FATAL_ERROR "A direct particle filter construction failed for another reason than the passkey.\n${output}")
    endif ()

    # A composed policy or a setter that does nothing is invisible at run time, so each anti-windup strategy
    # and set_integral must fail to compile with velocity_form, and for the stated reason.
    ctrlpp_try_canary(CTRLPP_PID_VELOCITY_BASELINE_COMPILED pid_velocity_anti_windup_canary.cpp)
    if (NOT CTRLPP_PID_VELOCITY_BASELINE_COMPILED)
        message(FATAL_ERROR
            "The pid velocity anti-windup canary is not a valid positive control: "
            "its velocity_form arm does not compile.\n"
            "${CTRLPP_PID_VELOCITY_BASELINE_COMPILED_OUTPUT}")
    endif ()
    foreach (strategy back_calc clamping conditional_integration)
        ctrlpp_try_canary(CTRLPP_PID_VELOCITY_${strategy}_COMPILED pid_velocity_anti_windup_canary.cpp
            -DCTRLPP_CANARY_VELOCITY_ANTI_WINDUP=ctrlpp::${strategy})
        set(output "${CTRLPP_PID_VELOCITY_${strategy}_COMPILED_OUTPUT}")
        if (CTRLPP_PID_VELOCITY_${strategy}_COMPILED)
            message(FATAL_ERROR "velocity_form with anti_windup<${strategy}> compiled.\n${output}")
        elseif (NOT output MATCHES "already prevents windup")
            message(FATAL_ERROR "velocity_form with anti_windup<${strategy}> failed for another reason.\n${output}")
        endif ()
    endforeach ()
    ctrlpp_try_canary(CTRLPP_PID_VELOCITY_SET_INTEGRAL_COMPILED pid_velocity_anti_windup_canary.cpp -DCTRLPP_CANARY_VELOCITY_SET_INTEGRAL)
    set(output "${CTRLPP_PID_VELOCITY_SET_INTEGRAL_COMPILED_OUTPUT}")
    if (CTRLPP_PID_VELOCITY_SET_INTEGRAL_COMPILED OR NOT output MATCHES "set_integral" OR NOT output MATCHES "constraint")
        message(FATAL_ERROR "set_integral on a velocity_form pid compiled, or failed for another reason.\n${output}")
    endif ()

    if (CTRLPP_CMAKE_FETCH_DEPS)
        FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.16.1
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
    # Catch2 shares the test TUs' -fno-rtti: a mixed program keeps whichever COMDAT
    # copy of an inline vtable it links first, and a -fno-rtti copy's null typeinfo
    # slot is then read by RTTI-enabled Catch2 code (UBSan: "invalid vptr").
    if (CTRLPP_CMAKE_FETCH_DEPS AND NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
        target_compile_definitions(Catch2 PUBLIC CATCH_CONFIG_DISABLE_EXCEPTIONS)
        target_compile_options(Catch2 PRIVATE -fno-exceptions -fno-rtti)
        target_compile_options(Catch2WithMain PRIVATE -fno-exceptions -fno-rtti)
    endif ()
    # RapidCheck needs RTTI for typeid and exceptions to report a falsified
    # property, so its tests can only link the exceptions carve-out's Catch2.
    if (CTRLPP_BUILD_PROPERTY_TESTS AND NOT CTRLPP_TESTS_WITH_EXCEPTIONS)
        message(FATAL_ERROR
            "Property tests (CTRLPP_BUILD_PROPERTY_TESTS=ON) need RTTI and exceptions, "
            "which the default -fno-exceptions -fno-rtti test tree does not provide. "
            "Configure them in the exceptions carve-out with -DCTRLPP_TESTS_WITH_EXCEPTIONS=ON.")
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
    # The source purity gate and the self-test that proves it can fail.
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
    add_subdirectory(integration)
endif ()

if (CTRLPP_BUILD_FUZZ_TESTS)
    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR "Fuzz targets require Clang with libFuzzer support")
    endif ()
    add_subdirectory(fuzz)
endif ()
