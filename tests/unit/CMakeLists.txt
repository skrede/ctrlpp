# Empty on MSVC, which needs /EHsc for its own runtime.
if (NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
    # Under -fno-exceptions Eigen 3.4's size-overflow check falls through to operator new(SIZE_MAX) instead of throwing,
    # so gcc keeps an impossible negative SimplicialLDLT workspace size live into aligned_malloc and reports it.
    set(CTRLPP_NOEXCEPT_FLAGS -fno-exceptions -fno-rtti $<$<CXX_COMPILER_ID:GNU>:-Wno-alloc-size-larger-than>)
else ()
    set(CTRLPP_NOEXCEPT_FLAGS "")
endif ()

function(make_test TEST_FILE)
    add_executable(${TEST_FILE} "${CMAKE_CURRENT_LIST_DIR}/${TEST_FILE}.cpp")
    target_link_libraries(${TEST_FILE} PRIVATE ctrlpp::ctrlpp Catch2::Catch2WithMain)
    target_compile_options(${TEST_FILE} PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME ${TEST_FILE} COMMAND ${TEST_FILE})
endfunction()

# For a test that needs exceptions: a REQUIRE_THROWS* assertion, or a link against a throwing solver.
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
make_test(pid_config_test)
make_test(pid_performance_test)
make_test(pid_mimo_test)
make_test(pid_output_test)
make_test(state_space_test)
make_test(observer_concept_test)
make_test(conversion_test)
make_test(analysis_test)
make_test(discretize_test)
make_test(pid_eigen_test)
make_test(dare_test)
make_test(care_test)
make_test(schur_reorder_test)
make_test(riccati_solution_test)
make_test(riccati_magnitude_test)
make_test(riccati_weight_scale_test)
make_test(psd_acceptance_test)
make_test(riccati_rank_deficient_test)
make_test(analysis_extended_test)
make_test(kalman_test)
make_test(luenberger_test)
make_test(lqr_test)
make_test(mpc_test)
make_test(output_admissible_set_test)
# The timeout is a runtime budget, not a tolerance: the explicit cases run on the order of the O_inf
# unit test's ~460 s in the Debug sanitizer tree, and the sweeps stay out of Debug trees.
add_executable(output_admissible_set_property_test "${CMAKE_CURRENT_LIST_DIR}/output_admissible_set_property_test.cpp")
target_link_libraries(output_admissible_set_property_test PRIVATE ctrlpp::ctrlpp Catch2::Catch2WithMain)
target_compile_options(output_admissible_set_property_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
add_test(NAME output_admissible_set_property_test COMMAND output_admissible_set_property_test "~[sweep]")
set_tests_properties(output_admissible_set_property_test PROPERTIES TIMEOUT 600)
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_test(NAME output_admissible_set_sweep_test COMMAND output_admissible_set_property_test "[sweep]")
    set_tests_properties(output_admissible_set_sweep_test PROPERTIES TIMEOUT 120)
endif ()
make_test(nmpc_test)
make_test(mpc_horizon_guard_test)
make_test(nmpc_static_solver_placement_test)
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
make_test(batch_arx_identity_test)
make_test(batch_arx_realization_test)
make_test(batch_arx_diagnostics_test)
make_test(moesp_test)
make_test(so3_test)
make_test(biquad_test)
make_test(fir_test)
make_test(complementary_filter_test)
make_test(mekf_test)
make_test(manifold_ukf_test)
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
make_test(cubic_spline_test)
make_test(online_planner_2nd_test)
make_test(online_planner_3rd_test)
make_test(online_planner_3rd_retarget_test)
# Checks the fuzz harness's oracle, included as "fuzz/online_planner_3rd_case.h".
target_include_directories(online_planner_3rd_retarget_test PRIVATE "${CMAKE_CURRENT_LIST_DIR}/..")
make_test(bspline_trajectory_test)
make_test(smoothing_spline_test)
make_test(synchronize_test)
make_test(mrac_core_test)
make_test(mrac_robustification_test)
make_test(mrac_mimo_test)
make_test(vector_biquad_test)
make_test_exceptions(l1_test)
make_test_exceptions(l1_projection_default_test)
make_test(pid_hardening_test)
make_test(lqr_hardening_test)
target_compile_definitions(lqr_hardening_test PRIVATE _GLIBCXX_ASSERTIONS)
make_test(dare_hardening_test)
make_test(dare_error_test)
make_test(dare_float_precision_test)
make_test(care_error_test)
make_test(care_convergence_anchor_test)
make_test(embedded_golden_tolerance_test)
make_test(embedded_telemetry_format_test)
# The reference's `quad` is `__float128`, which Apple clang on arm64 and MSVC lack; `long double` is
# binary64 there, so a substitute would stop being an oracle while still reporting verdicts. The probe
# compiles the type rather than trusting the macro alone.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
#if !defined(__SIZEOF_FLOAT128__)
#error no __float128
#endif
__float128 x = 0;
int main() { return static_cast<int>(x); }
" CTRLPP_HAVE_FLOAT128)

if (CTRLPP_HAVE_FLOAT128)
    make_test(care_quad_reference_test)
    # The binary128 reference lives beside the fuzz targets, included as "fuzz/care_quad_reference.h".
    target_include_directories(care_quad_reference_test PRIVATE "${CMAKE_CURRENT_LIST_DIR}/..")
    make_test(care_entitlement_test)
    target_include_directories(care_entitlement_test PRIVATE "${CMAKE_CURRENT_LIST_DIR}/..")
else ()
    message(STATUS "ctrlpp: no __float128 on this toolchain, skipping care_quad_reference_test and care_entitlement_test (the binary128 reference has no lower-precision substitute)")
endif ()
# MemorySanitizer's own strong operator new/delete collide at link time with nomalloc_harness.cpp's
# replacements (AddressSanitizer intercepts instead). The targets are skipped there because without the
# replacements the counter would see nothing and the allocation guarantee would pass untested.
string(FIND "${CMAKE_CXX_FLAGS}" "-fsanitize=memory" CTRLPP_MSAN_FLAG_POS)
if(CTRLPP_MSAN_FLAG_POS EQUAL -1)
    set(CTRLPP_NOMALLOC_TESTS_SUPPORTED TRUE)
else()
    set(CTRLPP_NOMALLOC_TESTS_SUPPORTED FALSE)
    message(STATUS "ctrlpp: MemorySanitizer detected, skipping *_nomalloc_test targets (operator new/delete replacement collides with libclang_rt.msan_cxx)")
endif()

if(CTRLPP_NOMALLOC_TESTS_SUPPORTED)
    add_library(ctrlpp_nomalloc_harness OBJECT "${CMAKE_CURRENT_LIST_DIR}/nomalloc_harness.cpp")
    target_link_libraries(ctrlpp_nomalloc_harness PRIVATE ctrlpp::ctrlpp)
    target_compile_options(ctrlpp_nomalloc_harness PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})

    # Every consumer of nomalloc_harness.h links its replacement allocation
    # functions exactly once and shares their process-global counter, so it runs
    # serially.
    function(make_nomalloc_test TEST_NAME)
        add_executable(${TEST_NAME} ${ARGN})
        target_link_libraries(${TEST_NAME} PRIVATE ctrlpp_nomalloc_harness ctrlpp::ctrlpp Catch2::Catch2WithMain)
        target_compile_options(${TEST_NAME} PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
        set_tests_properties(${TEST_NAME} PROPERTIES RUN_SERIAL TRUE)
    endfunction()

    # One estimator row's executable: its own unit plus PARTS units generated from
    # estimation_nomalloc_grid_part.cpp.in, each arming one cut of the row's
    # allocation-free grid. The cuts cover the grid exactly by construction, and
    # arm_part refuses a cut above max_points_per_unit at compile time.
    function(make_nomalloc_grid_test GRID_ROW GRID_LABEL GRID_TAG GRID_PARTS)
        set(test_name estimation_${GRID_ROW}_nomalloc_test)
        set(sources "${CMAKE_CURRENT_LIST_DIR}/${test_name}.cpp")
        math(EXPR last_part "${GRID_PARTS} - 1")
        foreach (GRID_PART RANGE ${last_part})
            math(EXPR GRID_PART_NUMBER "${GRID_PART} + 1")
            set(part_source "${CMAKE_CURRENT_BINARY_DIR}/${test_name}_part_${GRID_PART_NUMBER}.cpp")
            configure_file("${CMAKE_CURRENT_LIST_DIR}/estimation_nomalloc_grid_part.cpp.in" "${part_source}" @ONLY)
            list(APPEND sources "${part_source}")
        endforeach ()
        make_nomalloc_test(${test_name} ${sources})
        target_include_directories(${test_name} PRIVATE "${CMAKE_CURRENT_LIST_DIR}")
    endfunction()

    make_nomalloc_test(dare_care_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/dare_care_nomalloc_test.cpp")
    make_nomalloc_test(pid_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/pid_nomalloc_test.cpp")
    make_nomalloc_test(dsp_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/dsp_nomalloc_test.cpp")
    make_nomalloc_test(sysid_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/sysid_nomalloc_test.cpp")
    make_nomalloc_test(estimation_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/estimation_nomalloc_test.cpp")
    make_nomalloc_test(trajectory_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/trajectory_nomalloc_test.cpp")
    # A grid point costs its own instantiation chain, so each row's grid is cut into units of at most
    # max_points_per_unit points (estimation_nomalloc_grid.h); compiled as one unit the kalman_filter
    # row peaked at 12.3 GiB.
    make_nomalloc_grid_test(kalman kalman_filter kalman 5)
    make_nomalloc_grid_test(ekf ekf ekf 5)
    make_nomalloc_grid_test(ukf ukf ukf 6)
    make_nomalloc_grid_test(manifold_ukf manifold_ukf manifold_ukf 3)
    make_nomalloc_grid_test(mekf mekf mekf 4)
    # At the top of both axes these grids' whole-chain stack cost (docs/rt-safety-matrix.md) exceeds the
    # 1 MB main-thread stack MSVC reserves by default; the reserve is address space, committed on demand.
    if (MSVC)
        foreach (nomalloc_grid_target estimation_kalman_nomalloc_test estimation_ekf_nomalloc_test
                                      estimation_ukf_nomalloc_test estimation_manifold_ukf_nomalloc_test
                                      estimation_mekf_nomalloc_test)
            target_link_options(${nomalloc_grid_target} PRIVATE "/STACK:16777216")
        endforeach ()
    endif ()
endif()
make_test(lqr_continuous_test)
make_test(mrac_hardening_test)
make_test_exceptions(l1_hardening_test)
make_test(drift_admission_test)
make_test(one_output_gain_test)
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
make_test(trajectory_hardening_test)
make_test(so3_hardening_test)
make_test(state_space_hardening_test)
# The [!shouldfail] tag inverts a failing assertion into a pass, which relies on
# Catch2 catching the failure; under CATCH_CONFIG_DISABLE_EXCEPTIONS the failure
# path aborts instead, so this meta-test lives in the exceptions carve-out tree.
make_test_exceptions(shouldfail_smoke_test)
make_test(ukf_kf_equivalence_test)
make_test(anisotropic_p_anchor_test)
make_test(nees_nis_consistency_test)
make_test(dsp_response_anchor_test)
make_test(trajectory_envelope_anchor_test)
make_test(trajectory_rescale_anchor_test)
make_test(float_tier_test)
make_test(float_instantiation_test)
target_compile_options(float_tier_test PRIVATE ${CTRLPP_FLOAT_TIER_STRICT_FLAGS})
target_compile_options(float_instantiation_test PRIVATE ${CTRLPP_FLOAT_TIER_STRICT_FLAGS})
make_test(expected_test)
# The owned expected's throwing value() is asserted with CHECK_THROWS; its throw-free surface runs
# in the default tree through expected_test.
make_test_exceptions(expected_throws_test)
# The assignment matrix arms fixture members to throw from their copy and move constructors; CI builds
# the exceptions tree under the address and undefined-behavior sanitizers, so a double destruction fails.
make_test_exceptions(expected_transition_test)

# OSQP links a throwing solver.
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
    if(CTRLPP_NOMALLOC_TESTS_SUPPORTED)
        make_nomalloc_test(osqp_bounds_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/osqp_bounds_nomalloc_test.cpp")
        target_link_libraries(osqp_bounds_nomalloc_test PRIVATE ctrlpp::osqp)
    endif()
endif ()

# NLopt links a throwing solver.
if (CTRLPP_BUILD_NLOPT AND CTRLPP_TESTS_WITH_EXCEPTIONS)
    add_executable(nmpc_nlopt_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_nlopt_test.cpp")
    target_link_libraries(nmpc_nlopt_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain)
    target_compile_options(nmpc_nlopt_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nmpc_nlopt_test COMMAND nmpc_nlopt_test)

    find_package(Threads REQUIRED)
    add_executable(nlopt_containment_test "${CMAKE_CURRENT_LIST_DIR}/nlopt_containment_test.cpp")
    target_link_libraries(nlopt_containment_test PRIVATE ctrlpp::nlopt Catch2::Catch2WithMain Threads::Threads)
    target_compile_options(nlopt_containment_test PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME nlopt_containment_test COMMAND nlopt_containment_test)

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
    # nmpc_argmin_test and its AddressSanitizer twin assert that the throwing wrappers throw
    # (REQUIRE_THROWS*); the other argmin tests are throw-free.
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

        # Turns a dangling argmin back-pointer in the move-then-solve and copy-then-solve cases into a
        # hard failure.
        add_executable(nmpc_argmin_asan_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_argmin_test.cpp")
        target_link_libraries(nmpc_argmin_asan_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
        if (CTRLPP_BUILD_NLOPT)
            target_link_libraries(nmpc_argmin_asan_test PRIVATE ctrlpp::nlopt)
        endif ()
        if (CTRLPP_BUILD_OSQP)
            target_link_libraries(nmpc_argmin_asan_test PRIVATE ctrlpp::osqp)
            target_compile_definitions(nmpc_argmin_asan_test PRIVATE CTRLPP_HAS_OSQP)
        endif ()
        # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=116026, independent of Eigen, so this outlives the shared Eigen flags.
        target_compile_options(nmpc_argmin_asan_test PRIVATE ${CTRLPP_WARNING_FLAGS} -fsanitize=address -g
                               $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>)
        target_link_options(nmpc_argmin_asan_test PRIVATE -fsanitize=address)
        add_test(NAME nmpc_argmin_asan_test COMMAND nmpc_argmin_asan_test)
    endif ()

    add_executable(argmin_stepper_test "${CMAKE_CURRENT_LIST_DIR}/argmin_stepper_test.cpp")
    target_link_libraries(argmin_stepper_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_stepper_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_stepper_test COMMAND argmin_stepper_test)

    add_executable(argmin_budget_test "${CMAKE_CURRENT_LIST_DIR}/argmin_budget_test.cpp")
    target_link_libraries(argmin_budget_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_budget_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_budget_test COMMAND argmin_budget_test)

    add_executable(argmin_inequality_test "${CMAKE_CURRENT_LIST_DIR}/argmin_inequality_test.cpp")
    target_link_libraries(argmin_inequality_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_inequality_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_inequality_test COMMAND argmin_inequality_test)

    add_executable(argmin_setup_test "${CMAKE_CURRENT_LIST_DIR}/argmin_setup_test.cpp")
    target_link_libraries(argmin_setup_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_setup_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_setup_test COMMAND argmin_setup_test)

    add_executable(argmin_stop_criterion_test "${CMAKE_CURRENT_LIST_DIR}/argmin_stop_criterion_test.cpp")
    target_link_libraries(argmin_stop_criterion_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_stop_criterion_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_stop_criterion_test COMMAND argmin_stop_criterion_test)

    add_executable(nmhe_argmin_test "${CMAKE_CURRENT_LIST_DIR}/nmhe_argmin_test.cpp")
    target_link_libraries(nmhe_argmin_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(nmhe_argmin_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME nmhe_argmin_test COMMAND nmhe_argmin_test)

    # argmin_qp_solver.h defines CTRLPP_HAS_ARGMIN_QP when argmin's QP header is present, as it is in every argmin at
    # or above CTRLPP_ARGMIN_MIN_VERSION. Both backend tests gate their cases on it and keep one ungated case, since
    # Catch2 fails a binary that registers no test.
    add_executable(argmin_qp_solver_test "${CMAKE_CURRENT_LIST_DIR}/argmin_qp_solver_test.cpp")
    target_link_libraries(argmin_qp_solver_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    target_compile_options(argmin_qp_solver_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME argmin_qp_solver_test COMMAND argmin_qp_solver_test)

    add_executable(mpc_mhe_argmin_qp_test "${CMAKE_CURRENT_LIST_DIR}/mpc_mhe_argmin_qp_test.cpp")
    target_link_libraries(mpc_mhe_argmin_qp_test PRIVATE ctrlpp::argmin Catch2::Catch2WithMain)
    # The OSQP agreement legs need the throwing OSQP backend, so they join only where exceptions are on.
    if (CTRLPP_BUILD_OSQP AND CTRLPP_TESTS_WITH_EXCEPTIONS)
        target_link_libraries(mpc_mhe_argmin_qp_test PRIVATE ctrlpp::osqp)
        target_compile_definitions(mpc_mhe_argmin_qp_test PRIVATE CTRLPP_HAS_OSQP)
    endif ()
    target_compile_options(mpc_mhe_argmin_qp_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME mpc_mhe_argmin_qp_test COMMAND mpc_mhe_argmin_qp_test)

    if(CTRLPP_NOMALLOC_TESTS_SUPPORTED)
        make_nomalloc_test(nmpc_static_nomalloc_test "${CMAKE_CURRENT_LIST_DIR}/nmpc_static_nomalloc_test.cpp")
        target_link_libraries(nmpc_static_nomalloc_test PRIVATE ctrlpp::argmin)
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
