cmake_minimum_required(VERSION 3.28)

# Deliberately not PROJECT_IS_TOP_LEVEL: that is an ordinary variable set by
# project() and inherited by subdirectories, so this directory added by the
# repository root would read the root's own TRUE and take the standalone branch
# in both modes. Comparing the two source directories is the condition
# PROJECT_IS_TOP_LEVEL is itself derived from, and the only one that is correct
# before a project() call has run here.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(ctrlpp_validation)

    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    # Adding the repository root rather than the library directory: the library
    # listfile reads variables and include modules that only the root brings into
    # scope. OFF stops the nested root from adding this directory a second time.
    set(CTRLPP_BUILD_VALIDATION OFF)
    add_subdirectory(.. "${CMAKE_CURRENT_BINARY_DIR}/ctrlpp")

    # A variable set in a child directory does not travel back up to its parent,
    # so the root's warning set is read out of the scope that defines it.
    get_directory_property(CTRLPP_WARNING_FLAGS
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
        DEFINITION CTRLPP_WARNING_FLAGS
    )

    # CTest rather than a bare enable_testing(): the test runner reads a
    # configuration file only this module writes, and without it every run that
    # produces a results file opens by reporting that file as missing.
    include(CTest)
endif ()

# The cases compile under the same warning set as the rest of the tree rather
# than a weaker one of their own. The cost of that parity was measured before it
# was adopted: the case sources emit no diagnostics of their own under the full
# set on gcc and clang, and an optimized build under it emits no more than the
# three flags this harness used to carry. MSVC was not measured.
#
# That set promotes the discarded-result diagnostic to an error, and the reason
# weighs more here than anywhere else in the tree: these executables have no test
# framework and no assertions, so a dropped typed failure produces a
# plausible-looking CSV for a reference implementation to be compared against
# rather than a diagnostic anyone reads.
#
# Applied by the registration helper so every case inherits it, rather than
# repeated in eighteen leaf files where one could quietly drift.
function(ctrlpp_validation_case name)
    add_subdirectory(cases/${name})
    target_compile_options(${name} PRIVATE ${CTRLPP_WARNING_FLAGS})

    add_test(
        NAME ${name}
        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_case.sh" ${name} "$<TARGET_FILE:${name}>"
    )

    # Names the executable as a precondition of the test rather than letting the
    # wrapper discover it missing, so a case that was never built is reported as
    # a case that did not run and not as one that ran and disagreed.
    set_tests_properties(${name} PROPERTIES REQUIRED_FILES "$<TARGET_FILE:${name}>")
endfunction()

# Control module
ctrlpp_validation_case(dare_solution)
ctrlpp_validation_case(care_solution)
ctrlpp_validation_case(lqr_step_response)
ctrlpp_validation_case(lqr_finite_horizon)
ctrlpp_validation_case(lqi_tracking)
ctrlpp_validation_case(pid_linear_step)
ctrlpp_validation_case(pole_placement)

# Model module
ctrlpp_validation_case(c2d_zoh)
ctrlpp_validation_case(analysis_poles)
ctrlpp_validation_case(tf_ss_conversion)

# Estimation module
ctrlpp_validation_case(kalman_tracking)
ctrlpp_validation_case(luenberger_observer)

# DSP module
ctrlpp_validation_case(butterworth_filter)
ctrlpp_validation_case(fir_filter)

# Trajectory module
ctrlpp_validation_case(cubic_spline_natural)

# Lie groups module
ctrlpp_validation_case(so3_exp_log)

# System identification module
ctrlpp_validation_case(batch_arx_ident)
ctrlpp_validation_case(moesp_ident)
