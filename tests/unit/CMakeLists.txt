function(make_test TEST_FILE)
    add_executable(${TEST_FILE} "${CMAKE_CURRENT_LIST_DIR}/${TEST_FILE}.cpp")
    target_link_libraries(${TEST_FILE} PRIVATE ctrlpp::ctrlpp Catch2::Catch2WithMain)
    target_compile_options(${TEST_FILE} PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME ${TEST_FILE} COMMAND ${TEST_FILE})
endfunction()

make_test(pid_test)
make_test(state_space_test)
make_test(observer_concept_test)
make_test(conversion_test)
make_test(analysis_test)
make_test(discretise_test)
make_test(pid_eigen_test)
make_test(dare_test)
make_test(analysis_extended_test)
make_test(kalman_test)
make_test(luenberger_test)
make_test(lqr_test)
make_test(mpc_test)
make_test(nmpc_test)
make_test(numerical_diff_test)
make_test(runge_kutta_test)
make_test(ekf_test)
make_test(ukf_test)
make_test(particle_filter_test)

if (CTRLPP_BUILD_OSQP)
    add_executable(mpc_osqp_test "${CMAKE_CURRENT_LIST_DIR}/mpc_osqp_test.cpp")
    target_link_libraries(mpc_osqp_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mpc_osqp_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_osqp_test COMMAND mpc_osqp_test)

    add_executable(terminal_set_test "${CMAKE_CURRENT_LIST_DIR}/terminal_set_test.cpp")
    target_link_libraries(terminal_set_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(terminal_set_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME terminal_set_test COMMAND terminal_set_test)
endif ()

if (CTRLPP_BUILD_NLOPT)
    add_executable(nmpc_nlopt_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_nlopt_test.cpp")
    target_link_libraries(nmpc_nlopt_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmpc_nlopt_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmpc_nlopt_test COMMAND nmpc_nlopt_test)

    add_executable(nmpc_constraint_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_constraint_test.cpp")
    target_link_libraries(nmpc_constraint_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmpc_constraint_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmpc_constraint_test COMMAND nmpc_constraint_test)
endif ()

if (CTRLPP_BUILD_PROPERTY_TESTS)
    function(make_property_test TEST_FILE)
        add_executable(${TEST_FILE} "${CMAKE_CURRENT_LIST_DIR}/${TEST_FILE}.cpp")
        target_link_libraries(${TEST_FILE} PRIVATE
            ctrlpp::ctrlpp
            Catch2::Catch2WithMain
            rapidcheck
            rapidcheck_catch
        )
        target_compile_options(${TEST_FILE} PRIVATE ${CTRLPP_WARNING_FLAGS})
        add_test(NAME ${TEST_FILE} COMMAND ${TEST_FILE})
    endfunction()

    make_property_test(pid_property_test)
    make_property_test(dare_property_test)
    make_property_test(kalman_property_test)
endif ()

if (CTRLPP_BUILD_OSQP AND CTRLPP_BUILD_PROPERTY_TESTS)
    add_executable(mpc_property_test "${CMAKE_CURRENT_LIST_DIR}/mpc_property_test.cpp")
    target_link_libraries(mpc_property_test PRIVATE
        ctrlpp::osqp Catch2::Catch2WithMain rapidcheck rapidcheck_catch)
    target_compile_options(mpc_property_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_property_test COMMAND mpc_property_test)
endif ()

if (CTRLPP_BUILD_OSQP AND CTRLPP_BUILD_NLOPT)
    add_executable(mpc_closedloop_test "${CMAKE_CURRENT_LIST_DIR}/mpc_closedloop_test.cpp")
    target_link_libraries(mpc_closedloop_test PRIVATE
        ctrlpp::osqp ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(mpc_closedloop_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_closedloop_test COMMAND mpc_closedloop_test)
endif ()
