#include "ctrlpp/traj/trapezoidal_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// -- Test 1: Full trapezoidal profile ----------------------------------------
TEST_CASE("Trapezoidal: full trapezoidal profile", "[traj][trapezoidal]")
{
    // q0=0, q1=10, v_max=5, a_max=10
    // T_a = v_max/a_max = 0.5s, T_v = h/v_max - T_a = 2-0.5 = 1.5s, T_d = 0.5s, T = 2.5s
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(traj.duration());

    REQUIRE_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    REQUIRE_THAT(traj.peak_velocity(), WithinRel(5.0, 1e-12));
    REQUIRE(traj.is_triangular() == false);
}

// -- Test 2: Triangular degenerate -------------------------------------------
TEST_CASE("Trapezoidal: triangular degenerate case", "[traj][trapezoidal]")
{
    // q0=0, q1=1, v_max=10, a_max=2 -> h*a_max = 2 < v_max^2 = 100
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 1.0, .v_max = 10.0, .a_max = 2.0});

    REQUIRE(traj.is_triangular() == true);
    REQUIRE(traj.peak_velocity() < 10.0);
    REQUIRE(traj.peak_velocity() > 0.0);

    auto const pT = traj.evaluate(traj.duration());
    REQUIRE_THAT(pT.position[0], WithinAbs(1.0, 1e-12));
}

// -- Test 3: Velocity constraint ---------------------------------------------
TEST_CASE("Trapezoidal: velocity constraint satisfied", "[traj][trapezoidal]")
{
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

    double const eps = 1e-10;
    for (int i = 0; i <= 1000; ++i) {
        double const t = traj.duration() * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.velocity[0]) <= 5.0 + eps);
    }
}

// -- Test 4: Acceleration constraint -----------------------------------------
TEST_CASE("Trapezoidal: acceleration constraint satisfied", "[traj][trapezoidal]")
{
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

    double const eps = 1e-10;
    for (int i = 0; i <= 1000; ++i) {
        double const t = traj.duration() * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= 10.0 + eps);
    }
}

// -- Test 5: Phase boundary continuity ---------------------------------------
TEST_CASE("Trapezoidal: phase boundary continuity", "[traj][trapezoidal]")
{
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

    auto const phases = traj.phase_durations();
    double const T_a = phases[0];
    double const T_v = phases[1];

    // Boundary at T_a (accel -> cruise)
    auto const left1 = traj.evaluate(T_a - 1e-14);
    auto const right1 = traj.evaluate(T_a + 1e-14);
    REQUIRE_THAT(left1.position[0], WithinAbs(right1.position[0], 1e-10));
    REQUIRE_THAT(left1.velocity[0], WithinAbs(right1.velocity[0], 1e-10));

    // Boundary at T_a + T_v (cruise -> decel)
    auto const left2 = traj.evaluate(T_a + T_v - 1e-14);
    auto const right2 = traj.evaluate(T_a + T_v + 1e-14);
    REQUIRE_THAT(left2.position[0], WithinAbs(right2.position[0], 1e-10));
    REQUIRE_THAT(left2.velocity[0], WithinAbs(right2.velocity[0], 1e-10));
}

// -- Test 6: Non-null boundary conditions ------------------------------------
TEST_CASE("Trapezoidal: non-null boundary conditions", "[traj][trapezoidal]")
{
    ctrlpp::trapezoidal_trajectory<double> traj(
        {.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0, .v0 = 1.0, .v1 = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(traj.duration());

    REQUIRE_THAT(p0.velocity[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(pT.velocity[0], WithinAbs(2.0, 1e-12));
    REQUIRE_THAT(pT.position[0], WithinAbs(20.0, 1e-10));
}

// -- Test 7: Infeasible case - a_max adjusted per eq 3.15 --------------------
TEST_CASE("Trapezoidal: infeasible BCs adjust acceleration", "[traj][trapezoidal]")
{
    // a_max*h should be < |v0^2 - v1^2|/2 to trigger infeasibility
    // h=1, a_max=1, v0=0, v1=3 -> a_max*h=1 < (9-0)/2=4.5 -> infeasible
    ctrlpp::trapezoidal_trajectory<double> traj(
        {.q0 = 0.0, .q1 = 1.0, .v_max = 5.0, .a_max = 1.0, .v0 = 0.0, .v1 = 3.0});

    auto const pT = traj.evaluate(traj.duration());
    REQUIRE_THAT(pT.position[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pT.velocity[0], WithinAbs(3.0, 1e-10));
    REQUIRE(traj.duration() > 0.0);
}

// -- Test 8: Negative displacement -------------------------------------------
TEST_CASE("Trapezoidal: negative displacement", "[traj][trapezoidal]")
{
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 10.0, .q1 = 0.0, .v_max = 5.0, .a_max = 10.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(traj.duration());

    REQUIRE_THAT(p0.position[0], WithinAbs(10.0, 1e-12));
    REQUIRE_THAT(pT.position[0], WithinAbs(0.0, 1e-12));

    // velocity should be negative
    auto const p_mid = traj.evaluate(traj.duration() / 2.0);
    REQUIRE(p_mid.velocity[0] < 0.0);
}

// -- Test 9: Introspection ---------------------------------------------------
TEST_CASE("Trapezoidal: introspection", "[traj][trapezoidal]")
{
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

    auto const phases = traj.phase_durations();
    REQUIRE(phases.size() == 3);

    double const sum = phases[0] + phases[1] + phases[2];
    REQUIRE_THAT(sum, WithinAbs(traj.duration(), 1e-14));
}

// -- Test 10: Satisfies trajectory_segment concept ---------------------------
TEST_CASE("Trapezoidal: satisfies trajectory_segment concept", "[traj][trapezoidal]")
{
    STATIC_REQUIRE(ctrlpp::trajectory_segment<ctrlpp::trapezoidal_trajectory<double>, double, 1>);
    STATIC_REQUIRE(ctrlpp::trajectory_segment<ctrlpp::trapezoidal_trajectory<float>, float, 1>);
}
