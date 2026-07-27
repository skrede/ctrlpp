#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

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

// -- Test 11: Zero-displacement speed change is a typed rejection ------------
TEST_CASE("Trapezoidal: zero displacement with unequal boundary speeds is rejected",
          "[traj][trapezoidal]")
{
    // Both ramps of this family run toward one cruise velocity at or above each
    // boundary velocity, so the profile sweeps at least what the transition
    // between them already sweeps. eq. (3.15)'s remedy for a shorter command is
    // to divide the required speed change by the commanded displacement, which
    // at zero displacement asks for an unbounded acceleration; multiplied back
    // by the vanishing phase durations it produced, that lands as NaN in the
    // cruise duration rather than as a large finite profile.
    ctrlpp::trapezoidal_trajectory<double>::config const cfg{
        .q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = 0.0,
    };

    auto const created = ctrlpp::trapezoidal_trajectory<double>::try_create(cfg);
    REQUIRE_FALSE(created.has_value());
    REQUIRE(created.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);

    // The mirror command -- the same speed change, entered rather than left --
    // is the same infeasibility and carries the same value.
    auto const mirrored = ctrlpp::trapezoidal_trajectory<double>::try_create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.0, .v1 = 0.5});
    REQUIRE_FALSE(mirrored.has_value());
    REQUIRE(mirrored.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);

    // The plain constructor keeps working and reports a standstill, not a
    // profile whose duration, phases, and evaluation are all NaN.
    ctrlpp::trapezoidal_trajectory<double> traj(cfg);
    REQUIRE(std::isfinite(traj.duration()));
    REQUIRE_THAT(traj.duration(), WithinAbs(0.0, 0.0));
    for (auto const phase : traj.phase_durations()) {
        REQUIRE(std::isfinite(phase));
        REQUIRE(phase >= 0.0);
    }
    auto const p0 = traj.evaluate(0.0);
    REQUIRE_THAT(p0.position[0], WithinAbs(2.0, 0.0));
    REQUIRE_THAT(p0.velocity[0], WithinAbs(0.0, 0.0));
    REQUIRE_THAT(p0.acceleration[0], WithinAbs(0.0, 0.0));

    // Retiming a command the family cannot realize is a typed rejection too,
    // never arithmetic performed on the stationary stand-in.
    auto const rescaled = traj.rescale_to(1.0);
    REQUIRE_FALSE(rescaled.has_value());
    REQUIRE(rescaled.error() == ctrlpp::trajectory_error::unreachable_duration);
    REQUIRE_THAT(traj.duration(), WithinAbs(0.0, 0.0));
}

// -- Test 12: The rejection is confined to the infeasible set ----------------
TEST_CASE("Trapezoidal: zero displacement at equal boundary speeds is realizable",
          "[traj][trapezoidal]")
{
    using traj_t = ctrlpp::trapezoidal_trajectory<double>;

    // Equal boundary velocities: no speed change to cover, so the standstill is
    // the profile rather than a stand-in for a rejected one.
    auto const held = traj_t::try_create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = 0.5});
    REQUIRE(held.has_value());
    REQUIRE_THAT(held.value().duration(), WithinAbs(0.0, 0.0));

    // Opposed boundary velocities of equal magnitude: the ramp between them
    // sweeps exactly zero ground, so a zero command is realized by it exactly.
    auto const reversed = traj_t::try_create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = -0.5});
    REQUIRE(reversed.has_value());
    REQUIRE_THAT(reversed.value().duration(), WithinRel(1.0, 1e-15));

    // Rest to rest over zero displacement.
    auto const at_rest = traj_t::try_create({.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0});
    REQUIRE(at_rest.has_value());

    // A displacement small enough to force eq. (3.15)'s raise, but not so small
    // that the raised acceleration leaves the representable range: served, not
    // rejected. The boundary between the two is where the quotient overflows,
    // which is a property of the scalar type and is not a chosen threshold.
    auto const raised = traj_t::try_create(
        {.q0 = 0.0, .q1 = 1e-300, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = 0.0});
    REQUIRE(raised.has_value());
    REQUIRE(raised.value().duration() > 0.0);
    REQUIRE(std::isfinite(raised.value().duration()));
}
