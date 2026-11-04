cmake_minimum_required(VERSION 3.28)

# Not PROJECT_IS_TOP_LEVEL: it is an ordinary variable, so under the repository root this directory
# would inherit the root's TRUE. Comparing the source directories is what it is derived from, and is
# correct before any project() call here.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(ctrlpp_validation)

    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    # The repository root rather than lib/, because the library listfile reads variables and modules
    # only the root brings into scope. OFF stops the nested root adding this directory again.
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

# A waiver file is both sourced by a shell and read by a regular expression, so the quotes around a
# value are optional and stripped here.
function(ctrlpp_validation_waiver_value path key out)
    file(STRINGS "${path}" matches REGEX "^${key}=")

    set(value "")
    if (matches)
        list(GET matches 0 value)
        string(REGEX REPLACE "^${key}=[ \t]*" "" value "${value}")
        string(STRIP "${value}" value)
        string(REGEX REPLACE "^\"(.*)\"$" "\\1" value "${value}")
    endif ()

    set(${out} "${value}" PARENT_SCOPE)
endfunction()

# A waived case is a disabled test: it does not run, but keeps its row and its verbatim reason on the
# published results page, which is why a waiver without a reason and a date stops the configure.
function(ctrlpp_validation_waive_if_requested name)
    set(waiver "${CMAKE_CURRENT_SOURCE_DIR}/cases/${name}/waiver.cfg")
    if (NOT EXISTS "${waiver}")
        return()
    endif ()

    ctrlpp_validation_waiver_value("${waiver}" WAIVER_REASON reason)
    ctrlpp_validation_waiver_value("${waiver}" WAIVER_DATE date)
    if (reason STREQUAL "" OR date STREQUAL "")
        message(FATAL_ERROR
            "validation case ${name}: waiver.cfg must set both WAIVER_REASON and "
            "WAIVER_DATE, and neither may be empty"
        )
    endif ()

    set_tests_properties(${name} PROPERTIES DISABLED TRUE)
    message(STATUS "validation case ${name} waived ${date}: ${reason}")
endfunction()

# The cases compile under the tree's own warning set, applied here so that no leaf file can drift
# to a weaker one.
function(ctrlpp_validation_case name)
    add_subdirectory(cases/${name})
    target_compile_options(${name} PRIVATE ${CTRLPP_WARNING_FLAGS})

    add_test(
        NAME ${name}
        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_case.sh" ${name} "$<TARGET_FILE:${name}>"
    )

    # Reports a case that was never built as not run, rather than as one that ran and disagreed.
    set_tests_properties(${name} PROPERTIES REQUIRED_FILES "$<TARGET_FILE:${name}>")

    ctrlpp_validation_waive_if_requested(${name})
endfunction()

ctrlpp_validation_case(dare_solution)
ctrlpp_validation_case(care_solution)
ctrlpp_validation_case(lqr_step_response)
ctrlpp_validation_case(lqr_finite_horizon)
ctrlpp_validation_case(lqi_tracking)
ctrlpp_validation_case(pid_linear_step)
ctrlpp_validation_case(pole_placement)
ctrlpp_validation_case(lqr_closed_loop_settling)

ctrlpp_validation_case(c2d_zoh)
ctrlpp_validation_case(analysis_poles)
ctrlpp_validation_case(tf_ss_conversion)

ctrlpp_validation_case(kalman_tracking)
ctrlpp_validation_case(luenberger_observer)

ctrlpp_validation_case(butterworth_filter)
ctrlpp_validation_case(fir_filter)

ctrlpp_validation_case(cubic_spline_natural)

ctrlpp_validation_case(so3_exp_log)

ctrlpp_validation_case(batch_arx_ident)
ctrlpp_validation_case(moesp_ident)
