# The default (no-exceptions) tree dogfoods the throw-free discipline: ctrlpp's
# own test targets compile -fno-exceptions -fno-rtti. Empty in the exceptions
# carve-out tree and on MSVC (which needs /EHsc for its own runtime).
if (NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
    # Under -fno-exceptions Eigen 3.4's size-overflow check falls through to operator new(SIZE_MAX) instead of throwing,
    # so gcc keeps an impossible negative SimplicialLDLT workspace size live into aligned_malloc and reports it.
    # Remove once the fetched Eigen's no-exceptions allocation path no longer leaves that range reachable.
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
make_test(analysis_extended_test)
make_test(kalman_test)
make_test(luenberger_test)
make_test(lqr_test)
make_test(mpc_test)
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
# The reference's `quad` is `__float128`, which Apple clang on arm64 rejects and
# MSVC does not have. Substituting `long double` would be worse than skipping the
# target: it is binary64 on both of those toolchains, so the reference would carry
# no more precision than the answer it judges and would stop being an oracle while
# still reporting verdicts. The probe asserts the same condition the header does,
# and compiles the type rather than trusting the macro alone.
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
    # The binary128 reference lives beside the fuzz targets that consume it, and is
    # included as "fuzz/care_quad_reference.h"; this is what resolves that quote.
    target_include_directories(care_quad_reference_test PRIVATE "${CMAKE_CURRENT_LIST_DIR}/..")
else ()
    message(STATUS "ctrlpp: no __float128 on this toolchain, skipping care_quad_reference_test (the binary128 reference has no lower-precision substitute)")
endif ()
# MemorySanitizer ships its own strong operator new/delete in
# libclang_rt.msan_cxx, which collides at link time with the replacements in
# nomalloc_harness.cpp ("multiple definition of operator new"). AddressSanitizer
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
    # The five estimator rows whose frames a caller's dimensions size carry a
    # no-allocation dimension grid. A grid point costs its own instantiation
    # chain, so each row is one executable whose grid is cut into units of at
    # most max_points_per_unit points (estimation_nomalloc_grid.h); in one unit
    # the kalman_filter row peaked at 12.3 GiB. Measured one unit at a time with
    # gcc 16.2, the heaviest grid unit is ukf's, at 5.2 GiB at full -g and 3.7
    # GiB at -g1, below dare_care_nomalloc_test at 6.4 and 4.6.
    make_nomalloc_grid_test(kalman kalman_filter kalman 5)
    make_nomalloc_grid_test(ekf ekf ekf 5)
    make_nomalloc_grid_test(ukf ukf ukf 6)
    make_nomalloc_grid_test(manifold_ukf manifold_ukf manifold_ukf 3)
    make_nomalloc_grid_test(mekf mekf mekf 4)
    # These five arm dimensions whose whole-chain stack cost is published in
    # docs/rt-safety-matrix.md, and at the top of both axes that cost is above
    # the 1 MB main-thread stack MSVC reserves by default: the kalman_filter row
    # measures 71,224 bytes at 32 states against 4 outputs and grows
    # quadratically from there, so the 128-state points cost roughly a megabyte
    # on their own and the 128 x 128 corner several times that. The reserve
    # below is address space, committed on demand. Other toolchains here get 8 MB
    # from the operating system and need nothing.
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
# The owned expected's throwing value() accessor is asserted with CHECK_THROWS,
# which is unavailable under CATCH_CONFIG_DISABLE_EXCEPTIONS, so this TU lives in
# the exceptions carve-out tree; the throw-free surface and std interop run in the
# default tree via expected_test.
make_test_exceptions(expected_throws_test)
# The cross-state assignment matrix arms fixture members to throw from their copy
# and move constructors, so it needs exceptions; the carve-out tree is also the
# one CI builds under the address and undefined-behavior sanitizers, which is what
# makes a surviving double destruction a hard failure.
make_test_exceptions(expected_transition_test)

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
        # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=116026, independent of Eigen, so this outlives the shared Eigen
        # flags; remove once the supported gcc releases carry the fix.
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
    target_compile_options(mpc_mhe_argmin_qp_test PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_NOEXCEPT_FLAGS})
    add_test(NAME mpc_mhe_argmin_qp_test COMMAND mpc_mhe_argmin_qp_test)

    # Proof that the static NMPC path allocates nothing in steady state:
    # nmpc_static's solve is allocation-free at argmin's fixed-N floor. Consumes
    # the process-global allocation counter from nomalloc_harness.h, so it must
    # run serially, and is skipped under MemorySanitizer for the same
    # link-collision reason as the other *_nomalloc_test targets above.
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
