# The default (no-exceptions) tree dogfoods the throw-free discipline: ctrlpp's
# own test targets compile -fno-exceptions -fno-rtti. Empty in the exceptions
# carve-out tree and on MSVC (which needs /EHsc for its own runtime).
if (NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
    set(CTRLPP_NOEXCEPT_FLAGS -fno-exceptions -fno-rtti)
else ()
    set(CTRLPP_NOEXCEPT_FLAGS "")
endif ()

function(make_test TEST_FILE)
    add_executable(${TEST_FILE} "${CMAKE_CURRENT_LIST_DIR}/${TEST_FILE}.cpp")
    target_link_libraries(${TEST_FILE} PRIVATE ctrlpp::ctrlpp Catch2::Catch2WithMain)
    target_compile_options(${TEST_FILE} PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME ${TEST_FILE} COMMAND ${TEST_FILE})
endfunction()

# Registers a carve-out test that genuinely needs exceptions (a throwing-wrapper
# assertion via REQUIRE_THROWS*, or a link against a throwing solver). It builds
# ONLY in the exceptions tree (CTRLPP_TESTS_WITH_EXCEPTIONS=ON), where Catch2
# keeps exceptions on; in the default -fno-exceptions tree it is a no-op.
function(make_test_exceptions TEST_FILE)
    if (NOT CTRLPP_TESTS_WITH_EXCEPTIONS)
        return()
    endif ()
    add_executable(${TEST_FILE} "${CMAKE_CURRENT_LIST_DIR}/${TEST_FILE}.cpp")
    target_link_libraries(${TEST_FILE} PRIVATE ctrlpp::ctrlpp Catch2::Catch2WithMain)
    target_compile_options(${TEST_FILE} PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME ${TEST_FILE} COMMAND ${TEST_FILE})
endfunction()

make_test(pid_core_test)
make_test(pid_integration_method_test)
make_test(pid_derivative_test)
make_test(pid_filtering_test)
make_test(pid_antiwindup_test)
make_test(pid_velocity_form_test)
make_test(pid_tracking_test)
make_test(pid_gain_scheduling_test)
make_test(pid_performance_test)
make_test(pid_mimo_test)
make_test(pid_output_test)
make_test(state_space_test)
make_test(observer_concept_test)
make_test(conversion_test)
make_test(analysis_test)
make_test(discretise_test)
make_test(pid_eigen_test)
make_test(dare_test)
make_test(care_test)
make_test(schur_reorder_test)
make_test(riccati_solution_test)
make_test(analysis_extended_test)
make_test(kalman_test)
make_test(luenberger_test)
make_test(lqr_test)
make_test(mpc_test)
make_test(nmpc_test)
make_test(numerical_diff_test)
make_test(ekf_convergence_test)
make_test(ekf_covariance_test)
make_test(ukf_sigma_test)
make_test(ukf_tracking_test)
make_test(ukf_strategy_test)
make_test(particle_filter_resampling_test)
make_test(particle_filter_tracking_test)
make_test(particle_filter_config_test)
make_test(fit_metrics_test)
make_test(rls_test)
make_test(recursive_arx_test)
make_test(batch_arx_test)
make_test(moesp_test)
make_test(so3_test)
make_test(biquad_test)
make_test(fir_test)
make_test_exceptions(complementary_filter_test)
make_test_exceptions(mekf_test)
make_test_exceptions(manifold_ukf_test)
make_test(trajectory_types_test)
make_test(path_test)
make_test(polynomial_trajectory_test)
make_test(trajectory_adapter_test)
make_test(piecewise_trajectory_test)
make_test(trapezoidal_trajectory_test)
make_test(double_s_trajectory_test)
make_test(modified_trap_trajectory_test)
make_test(modified_sin_trajectory_test)
make_test(time_scaling_test)
make_test_exceptions(cubic_spline_test)
make_test_exceptions(online_planner_2nd_test)
make_test_exceptions(online_planner_3rd_test)
make_test_exceptions(bspline_trajectory_test)
make_test_exceptions(smoothing_spline_test)
make_test(synchronize_test)
make_test(mrac_core_test)
make_test(mrac_robustification_test)
make_test(mrac_mimo_test)
make_test(vector_biquad_test)
make_test_exceptions(l1_test)
make_test_exceptions(l1_projection_default_test)
make_test(pid_hardening_test)
make_test(lqr_hardening_test)
make_test(dare_hardening_test)
make_test(dare_error_test)
make_test(care_error_test)
# MemorySanitizer ships its own strong operator new/delete in
# libclang_rt.msan_cxx, which collides at link time with the replacements in
# nomalloc_harness.h ("multiple definition of operator new"). AddressSanitizer
# intercepts allocation instead of replacing these operators, so only MSan is
# affected. These targets are skipped there rather than having the harness drop
# its replacements under MSan: without them the counter would see nothing and
# the allocation guarantee would silently pass instead of being tested.
string(FIND "${CMAKE_CXX_FLAGS}" "-fsanitize=memory" CTRLPP_MSAN_FLAG_POS)
if(CTRLPP_MSAN_FLAG_POS EQUAL -1)
    set(CTRLPP_NOMALLOC_TESTS_SUPPORTED TRUE)
else()
    set(CTRLPP_NOMALLOC_TESTS_SUPPORTED FALSE)
    message(STATUS "ctrlpp: MemorySanitizer detected, skipping *_nomalloc_test targets (operator new/delete replacement collides with libclang_rt.msan_cxx)")
endif()

if(CTRLPP_NOMALLOC_TESTS_SUPPORTED)
    make_test(dare_care_nomalloc_test)
    make_test(pid_nomalloc_test)
    make_test(dsp_nomalloc_test)
    make_test(sysid_nomalloc_test)
    make_test(estimation_nomalloc_test)
    make_test(trajectory_nomalloc_test)
    # Every *_nomalloc_test consumes the process-global allocation counter from
    # nomalloc_harness.h and must be registered RUN_SERIAL.
    set_tests_properties(dare_care_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(pid_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(dsp_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(sysid_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(estimation_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(trajectory_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
endif()
make_test(lqr_continuous_test)
make_test(mrac_hardening_test)
make_test_exceptions(l1_hardening_test)
make_test(kalman_hardening_test)
make_test(ekf_hardening_test)
make_test(ukf_hardening_test)
make_test(particle_filter_hardening_test)
make_test(luenberger_hardening_test)
make_test_exceptions(complementary_filter_hardening_test)
make_test_exceptions(mekf_hardening_test)
make_test_exceptions(manifold_ukf_hardening_test)
make_test(sysid_hardening_test)
make_test(dsp_hardening_test)
make_test_exceptions(trajectory_hardening_test)
make_test(so3_hardening_test)
make_test(state_space_hardening_test)
# The [!shouldfail] tag inverts a failing assertion into a pass, which relies on
# Catch2 catching the failure; under CATCH_CONFIG_DISABLE_EXCEPTIONS the failure
# path aborts instead, so this meta-test lives in the exceptions carve-out tree.
make_test_exceptions(shouldfail_smoke_test)
make_test(ukf_kf_equivalence_test)
make_test_exceptions(anisotropic_p_anchor_test)
make_test_exceptions(nees_nis_consistency_test)
make_test(dsp_response_anchor_test)
make_test(trajectory_envelope_anchor_test)
make_test(float_tier_test)
make_test_exceptions(float_instantiation_test)
make_test(expected_test)
# The owned expected's throwing value() accessor is asserted with CHECK_THROWS,
# which is unavailable under CATCH_CONFIG_DISABLE_EXCEPTIONS, so this TU lives in
# the exceptions carve-out tree; the throw-free surface and std interop run in the
# default tree via expected_test.
make_test_exceptions(expected_throws_test)

# OSQP links a throwing solver, so these TUs belong to the exceptions carve-out
# tree and cannot compile in the default -fno-exceptions tree.
if (CTRLPP_BUILD_OSQP AND CTRLPP_TESTS_WITH_EXCEPTIONS)
    add_executable(mpc_osqp_test "${CMAKE_CURRENT_LIST_DIR}/mpc_osqp_test.cpp")
    target_link_libraries(mpc_osqp_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mpc_osqp_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_osqp_test COMMAND mpc_osqp_test)

    add_executable(terminal_set_test "${CMAKE_CURRENT_LIST_DIR}/terminal_set_test.cpp")
    target_link_libraries(terminal_set_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(terminal_set_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME terminal_set_test COMMAND terminal_set_test)

    add_executable(mhe_test "${CMAKE_CURRENT_LIST_DIR}/mhe_test.cpp")
    target_link_libraries(mhe_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mhe_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mhe_test COMMAND mhe_test)

    add_executable(mhe_kf_equivalence_test "${CMAKE_CURRENT_LIST_DIR}/mhe_kf_equivalence_test.cpp")
    target_link_libraries(mhe_kf_equivalence_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mhe_kf_equivalence_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mhe_kf_equivalence_test COMMAND mhe_kf_equivalence_test)

    add_executable(mpc_hardening_test "${CMAKE_CURRENT_LIST_DIR}/mpc_hardening_test.cpp")
    target_link_libraries(mpc_hardening_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mpc_hardening_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_hardening_test COMMAND mpc_hardening_test)

    add_executable(mhe_hardening_test "${CMAKE_CURRENT_LIST_DIR}/mhe_hardening_test.cpp")
    target_link_libraries(mhe_hardening_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mhe_hardening_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mhe_hardening_test COMMAND mhe_hardening_test)

    add_executable(mpc_coverage_test "${CMAKE_CURRENT_LIST_DIR}/mpc_coverage_test.cpp")
    target_link_libraries(mpc_coverage_test PRIVATE ctrlpp::osqp Catch2::Catch2WithMain)
    target_compile_options(mpc_coverage_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_coverage_test COMMAND mpc_coverage_test)
endif ()

# NLopt links a throwing solver, so these TUs belong to the exceptions carve-out
# tree and cannot compile in the default -fno-exceptions tree.
if (CTRLPP_BUILD_NLOPT AND CTRLPP_TESTS_WITH_EXCEPTIONS)
    add_executable(nmpc_nlopt_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_nlopt_test.cpp")
    target_link_libraries(nmpc_nlopt_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmpc_nlopt_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmpc_nlopt_test COMMAND nmpc_nlopt_test)

    add_executable(nmpc_constraint_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_constraint_test.cpp")
    target_link_libraries(nmpc_constraint_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmpc_constraint_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmpc_constraint_test COMMAND nmpc_constraint_test)

    add_executable(nmhe_test "${CMAKE_CURRENT_LIST_DIR}/nmhe_test.cpp")
    target_link_libraries(nmhe_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmhe_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmhe_test COMMAND nmhe_test)

    add_executable(nmpc_coverage_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_coverage_test.cpp")
    target_link_libraries(nmpc_coverage_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmpc_coverage_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmpc_coverage_test COMMAND nmpc_coverage_test)
endif ()

if (CTRLPP_BUILD_ARGMIN)
    # nmpc_argmin_test asserts the throwing wrappers throw (REQUIRE_THROWS*), so
    # it and its AddressSanitizer twin live in the exceptions carve-out tree; the
    # remaining argmin tests are throw-free and run in the default tree.
    if (CTRLPP_TESTS_WITH_EXCEPTIONS)
        add_executable(nmpc_argmin_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_argmin_test.cpp")
        target_link_libraries(nmpc_argmin_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
        if (CTRLPP_BUILD_NLOPT)
            target_link_libraries(nmpc_argmin_test PRIVATE ctrlpp::nlopt)
        endif ()
        if (CTRLPP_BUILD_OSQP)
            target_link_libraries(nmpc_argmin_test PRIVATE ctrlpp::osqp)
            target_compile_definitions(nmpc_argmin_test PRIVATE CTRLPP_HAS_OSQP)
        endif ()
        target_compile_options(nmpc_argmin_test PRIVATE ${CTRLPP_WARNING_FLAGS})
        add_test(NAME nmpc_argmin_test COMMAND nmpc_argmin_test)

        # AddressSanitizer variant of the move/copy safety suite. Same source, built
        # with -fsanitize=address so the move-then-solve / copy-then-solve cases turn
        # a dangling argmin back-pointer into a hard ASAN failure.
        add_executable(nmpc_argmin_asan_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_argmin_test.cpp")
        target_link_libraries(nmpc_argmin_asan_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
        if (CTRLPP_BUILD_NLOPT)
            target_link_libraries(nmpc_argmin_asan_test PRIVATE ctrlpp::nlopt)
        endif ()
        if (CTRLPP_BUILD_OSQP)
            target_link_libraries(nmpc_argmin_asan_test PRIVATE ctrlpp::osqp)
            target_compile_definitions(nmpc_argmin_asan_test PRIVATE CTRLPP_HAS_OSQP)
        endif ()
        target_compile_options(nmpc_argmin_asan_test PRIVATE ${CTRLPP_WARNING_FLAGS} -fsanitize=address -g)
        target_link_options(nmpc_argmin_asan_test PRIVATE -fsanitize=address)
        add_test(NAME nmpc_argmin_asan_test COMMAND nmpc_argmin_asan_test)
    endif ()

    add_executable(argmin_stepper_test "${CMAKE_CURRENT_LIST_DIR}/argmin_stepper_test.cpp")
    target_link_libraries(argmin_stepper_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_stepper_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_stepper_test COMMAND argmin_stepper_test)

    add_executable(nmhe_argmin_test "${CMAKE_CURRENT_LIST_DIR}/nmhe_argmin_test.cpp")
    target_link_libraries(nmhe_argmin_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(nmhe_argmin_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME nmhe_argmin_test COMMAND nmhe_argmin_test)

    # QP backend policy over argmin's sparse ADMM solver. Self-gates on the
    # argmin QP header via __has_include, so it is a trivially-passing empty TU
    # at the production pin (which predates argmin's QP headers) and exercises
    # the concept + a functional solve against a checkout that ships them.
    add_executable(argmin_qp_solver_test "${CMAKE_CURRENT_LIST_DIR}/argmin_qp_solver_test.cpp")
    target_link_libraries(argmin_qp_solver_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_qp_solver_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_qp_solver_test COMMAND argmin_qp_solver_test)

    # EMB-02 proof for the static NMPC path: nmpc_static's steady-state solve is
    # allocation-free at argmin's fixed-N floor. Consumes the process-global
    # allocation counter from nomalloc_harness.h, so it must run serially, and
    # is skipped under MemorySanitizer for the same link-collision reason as the
    # other *_nomalloc_test targets above.
    if(CTRLPP_NOMALLOC_TESTS_SUPPORTED)
        add_executable(nmpc_static_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_static_nomalloc_test.cpp")
        target_link_libraries(nmpc_static_nomalloc_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
        target_compile_options(nmpc_static_nomalloc_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
        add_test(NAME nmpc_static_nomalloc_test COMMAND nmpc_static_nomalloc_test)
        set_tests_properties(nmpc_static_nomalloc_test PROPERTIES RUN_SERIAL TRUE)
    endif()
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
    make_property_test(so3_property_test)
    make_property_test(estimator_property_test)
    make_property_test(dsp_property_test)
    make_property_test(sysid_property_test)
endif ()

if (CTRLPP_BUILD_OSQP AND CTRLPP_BUILD_PROPERTY_TESTS AND CTRLPP_TESTS_WITH_EXCEPTIONS)
    add_executable(mpc_property_test "${CMAKE_CURRENT_LIST_DIR}/mpc_property_test.cpp")
    target_link_libraries(mpc_property_test PRIVATE
        ctrlpp::osqp Catch2::Catch2WithMain rapidcheck rapidcheck_catch)
    target_compile_options(mpc_property_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_property_test COMMAND mpc_property_test)
endif ()

if (CTRLPP_BUILD_OSQP AND CTRLPP_BUILD_NLOPT AND CTRLPP_TESTS_WITH_EXCEPTIONS)
    add_executable(mpc_closedloop_test "${CMAKE_CURRENT_LIST_DIR}/mpc_closedloop_test.cpp")
    target_link_libraries(mpc_closedloop_test PRIVATE
        ctrlpp::osqp ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(mpc_closedloop_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME mpc_closedloop_test COMMAND mpc_closedloop_test)
endif ()
