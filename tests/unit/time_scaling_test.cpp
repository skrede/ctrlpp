#include "ctrlpp/trajectory/time_scaling.h"
#include "ctrlpp/trajectory/cubic_path.h"
#include "ctrlpp/trajectory/cycloidal_path.h"
#include "ctrlpp/trajectory/septic_path.h"
#include "ctrlpp/trajectory/trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// -- Test 1: Cycloidal path, velocity-limited --------------------------------
TEST_CASE("Time scaling: cycloidal path velocity-limited", "[traj][time_scaling]")
{
    // h=10, v_max=5, a_max=100, j_max=1000
    // peak_derivs = {2.0, 2*pi, 4*pi^2}
    // T_vel = 10*2/5 = 4.0
    // T_acc = sqrt(10*2*pi/100) = sqrt(0.6283) ~ 0.7927
    // T_jrk = cbrt(10*4*pi^2/1000) = cbrt(0.3948) ~ 0.7338
    // T = max(4.0, 0.7927, 0.7338) = 4.0

    auto const pd = ctrlpp::cycloidal_path_peak_derivatives<double>();
    auto const T = ctrlpp::compute_min_duration(10.0, pd, 5.0, 100.0, 1000.0);

    REQUIRE_THAT(T, WithinRel(4.0, 1e-10));
}

// -- Test 2: Cubic path, acceleration-limited --------------------------------
TEST_CASE("Time scaling: cubic path acceleration-limited", "[traj][time_scaling]")
{
    // h=1, v_max=100, a_max=2, j_max=1000
    // peak_derivs = {1.5, 6.0, 12.0}
    // T_vel = 1*1.5/100 = 0.015
    // T_acc = sqrt(1*6.0/2) = sqrt(3.0) ~ 1.7321
    // T_jrk = cbrt(1*12.0/1000) = cbrt(0.012) ~ 0.2289
    // T = T_acc

    auto const pd = ctrlpp::cubic_path_peak_derivatives<double>();
    auto const T = ctrlpp::compute_min_duration(1.0, pd, 100.0, 2.0, 1000.0);

    REQUIRE_THAT(T, WithinRel(std::sqrt(3.0), 1e-10));
}

// -- Test 3: Septic path, jerk-limited ---------------------------------------
TEST_CASE("Time scaling: septic path jerk-limited", "[traj][time_scaling]")
{
    // h=1, v_max=100, a_max=100, j_max=1
    // peak_derivs = {35/16, 7.5132, 52.5}
    // T_vel = 1*(35/16)/100 = 0.021875
    // T_acc = sqrt(1*7.5132/100) = sqrt(0.075132) ~ 0.2741
    // T_jrk = cbrt(1*52.5/1) = cbrt(52.5) ~ 3.7476
    // T = T_jrk

    auto const pd = ctrlpp::septic_path_peak_derivatives<double>();
    auto const T = ctrlpp::compute_min_duration(1.0, pd, 100.0, 100.0, 1.0);

    REQUIRE_THAT(T, WithinRel(std::cbrt(52.5), 1e-10));
}

// -- Test 4: Verify resulting trajectory satisfies limits ---------------------
TEST_CASE("Time scaling: trajectory respects limits", "[traj][time_scaling]")
{
    double const h = 5.0;
    double const v_max = 3.0;
    double const a_max = 10.0;
    double const j_max = 1000.0;

    auto const pd = ctrlpp::cycloidal_path_peak_derivatives<double>();
    auto const T = ctrlpp::compute_min_duration(h, pd, v_max, a_max, j_max);

    // Create trajectory with computed T
    Eigen::Matrix<double, 1, 1> const q0 = Eigen::Matrix<double, 1, 1>::Zero();
    Eigen::Matrix<double, 1, 1> const q1 = Eigen::Matrix<double, 1, 1>::Constant(h);
    auto const traj = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, T);

    double const eps = 1e-8;
    for (int i = 0; i <= 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.velocity[0]) <= v_max + eps);
        REQUIRE(std::abs(pt.acceleration[0]) <= a_max + eps);
    }
}

// -- Test 5: Vector overload - per-DOF durations -----------------------------
TEST_CASE("Time scaling: vector overload per-DOF", "[traj][time_scaling]")
{
    auto const pd = ctrlpp::cycloidal_path_peak_derivatives<double>();
    std::array<double, 3> const h = {1.0, 5.0, 10.0};

    auto const T_per_dof = ctrlpp::compute_min_duration(h, pd, 5.0, 100.0, 1000.0);

    // All three should be different (different h -> different T)
    REQUIRE(T_per_dof[0] != T_per_dof[1]);
    REQUIRE(T_per_dof[1] != T_per_dof[2]);

    // Each should equal the scalar overload
    for (std::size_t i = 0; i < 3; ++i) {
        auto const T_scalar = ctrlpp::compute_min_duration(h[i], pd, 5.0, 100.0, 1000.0);
        REQUIRE_THAT(T_per_dof[i], WithinAbs(T_scalar, 1e-14));
    }
}

// -- Test 6: Scalar synchronized overload ------------------------------------
TEST_CASE("Time scaling: scalar synchronized overload", "[traj][time_scaling]")
{
    auto const pd = ctrlpp::cycloidal_path_peak_derivatives<double>();
    std::array<double, 3> const h = {1.0, 5.0, 10.0};

    auto const T_sync = ctrlpp::compute_min_duration_sync(h, pd, 5.0, 100.0, 1000.0);
    auto const T_per_dof = ctrlpp::compute_min_duration(h, pd, 5.0, 100.0, 1000.0);

    // Sync should be the maximum of per-DOF
    auto const T_max = *std::max_element(T_per_dof.begin(), T_per_dof.end());
    REQUIRE_THAT(T_sync, WithinAbs(T_max, 1e-14));
}

// -- Test 7: Composite profiles do NOT use time_scaling (conceptual) ---------
TEST_CASE("Time scaling: composite profiles solve timing internally", "[traj][time_scaling]")
{
    // Per D-14: Composite profiles (trapezoidal, double-S) do NOT use time_scaling.
    // They solve timing internally from v_max, a_max constraints.
    // This test documents the design decision -- no code dependency to verify.
    //
    // time_scaling is for elementary paths lifted into trajectories via
    // make_trajectory(path, q0, q1, T), where T = compute_min_duration(...).

    SUCCEED("Conceptual test: time_scaling is for elementary paths only (D-14)");
}
