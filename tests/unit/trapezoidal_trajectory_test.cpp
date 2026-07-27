#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

// create() is the only construction path and it is fallible, so every profile a
// test uses is built through this helper, which asserts the command was
// realizable rather than letting a rejection pass as a profile. The rejection
// cases below never come through here: they assert the specific enumerator.
auto built(ctrlpp::trapezoidal_trajectory<double>::config const& cfg)
    -> ctrlpp::trapezoidal_trajectory<double>
{
    auto created = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

}

// -- Test 1: Full trapezoidal profile ----------------------------------------
TEST_CASE("Trapezoidal: full trapezoidal profile", "[traj][trapezoidal]")
{
    // q0=0, q1=10, v_max=5, a_max=10
    // T_a = v_max/a_max = 0.5s, T_v = h/v_max - T_a = 2-0.5 = 1.5s, T_d = 0.5s, T = 2.5s
    auto const traj = built({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

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
    auto const traj = built({.q0 = 0.0, .q1 = 1.0, .v_max = 10.0, .a_max = 2.0});

    REQUIRE(traj.is_triangular() == true);
    REQUIRE(traj.peak_velocity() < 10.0);
    REQUIRE(traj.peak_velocity() > 0.0);

    auto const pT = traj.evaluate(traj.duration());
    REQUIRE_THAT(pT.position[0], WithinAbs(1.0, 1e-12));
}

// -- Test 3: Velocity constraint ---------------------------------------------
TEST_CASE("Trapezoidal: velocity constraint satisfied", "[traj][trapezoidal]")
{
    auto const traj = built({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

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
    auto const traj = built({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

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
    auto const traj = built({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

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
    auto const traj = built({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0, .v0 = 1.0, .v1 = 2.0});

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
    auto const traj = built({.q0 = 0.0, .q1 = 1.0, .v_max = 5.0, .a_max = 1.0, .v0 = 0.0, .v1 = 3.0});

    auto const pT = traj.evaluate(traj.duration());
    REQUIRE_THAT(pT.position[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pT.velocity[0], WithinAbs(3.0, 1e-10));
    REQUIRE(traj.duration() > 0.0);
}

// -- Test 8: Negative displacement -------------------------------------------
TEST_CASE("Trapezoidal: negative displacement", "[traj][trapezoidal]")
{
    auto const traj = built({.q0 = 10.0, .q1 = 0.0, .v_max = 5.0, .a_max = 10.0});

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
    auto const traj = built({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

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
    // at zero displacement asks for an unbounded acceleration; every quantity
    // downstream is that value multiplied by a phase duration derived from it,
    // so it propagates as NaN rather than as a large finite profile. The
    // rejection is decided from the remedy's own representability, before
    // anything is built.
    ctrlpp::trapezoidal_trajectory<double>::config const cfg{
        .q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = 0.0,
    };

    auto const created = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
    REQUIRE_FALSE(created.has_value());
    REQUIRE(created.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);

    // The mirror command -- the same speed change, entered rather than left --
    // is the same infeasibility and carries the same value.
    auto const mirrored = ctrlpp::trapezoidal_trajectory<double>::create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.0, .v1 = 0.5});
    REQUIRE_FALSE(mirrored.has_value());
    REQUIRE(mirrored.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);
}

// -- Test 12: The rejection is confined to the infeasible set ----------------
TEST_CASE("Trapezoidal: zero displacement at equal boundary speeds is realizable",
          "[traj][trapezoidal]")
{
    using traj_t = ctrlpp::trapezoidal_trajectory<double>;

    // Equal boundary velocities: no speed change to cover, so the standstill is
    // the profile rather than a stand-in for a rejected one.
    auto const held = traj_t::create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = 0.5});
    REQUIRE(held.has_value());
    REQUIRE_THAT(held.value().duration(), WithinAbs(0.0, 0.0));

    // Opposed boundary velocities of equal magnitude: the ramp between them
    // sweeps exactly zero ground, so a zero command is realized by it exactly.
    auto const reversed = traj_t::create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = -0.5});
    REQUIRE(reversed.has_value());
    REQUIRE_THAT(reversed.value().duration(), WithinRel(1.0, 1e-15));

    // Rest to rest over zero displacement.
    auto const at_rest = traj_t::create({.q0 = 2.0, .q1 = 2.0, .v_max = 2.0, .a_max = 1.0});
    REQUIRE(at_rest.has_value());

    // A displacement small enough to force eq. (3.15)'s raise, but not so small
    // that the raised acceleration leaves the representable range: served, not
    // rejected. The boundary between the two is where the quotient overflows,
    // which is a property of the scalar type and is not a chosen threshold.
    auto const raised = traj_t::create(
        {.q0 = 0.0, .q1 = 1e-300, .v_max = 2.0, .a_max = 1.0, .v0 = 0.5, .v1 = 0.0});
    REQUIRE(raised.has_value());
    REQUIRE(raised.value().duration() > 0.0);
    REQUIRE(std::isfinite(raised.value().duration()));
}

// -- Test 13: the velocity limit is a precondition on the boundary velocities --
TEST_CASE("Trapezoidal: a boundary velocity above the velocity limit is rejected",
          "[traj][trapezoidal]")
{
    using traj_t = ctrlpp::trapezoidal_trajectory<double>;

    // The acceleration phase spans (v_v - v0) / a with the cruise velocity held
    // at or under the limit, so a boundary velocity above the limit makes that
    // span negative: the profile would have to run backwards in time. The
    // alternative -- raising the limit to whatever the caller passed as a
    // boundary velocity -- returns a profile that violates a bound the caller
    // stated, so the command is rejected instead.
    auto const entering = traj_t::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 1.0, .a_max = 1.0, .v0 = 5.0, .v1 = 0.0});
    REQUIRE_FALSE(entering.has_value());
    REQUIRE(entering.error() == ctrlpp::trajectory_error::boundary_velocity_exceeds_limit);

    // The limit bounds the SPEED, so the sign of the boundary velocity does not
    // enter and the terminal one is held to it identically.
    auto const leaving = traj_t::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 1.0, .a_max = 1.0, .v0 = 0.0, .v1 = -5.0});
    REQUIRE_FALSE(leaving.has_value());
    REQUIRE(leaving.error() == ctrlpp::trajectory_error::boundary_velocity_exceeds_limit);

    // Reaching the limit exactly is inside the domain, not outside it: the ramp
    // attached to that boundary velocity vanishes, which is a degenerate shape
    // rather than an impossible one.
    auto const at_limit = traj_t::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 1.0, .a_max = 1.0, .v0 = 1.0, .v1 = 0.0});
    REQUIRE(at_limit.has_value());
    auto const phases = at_limit.value().phase_durations();
    REQUIRE_THAT(phases[0], WithinAbs(0.0, 0.0));
    REQUIRE(at_limit.value().duration() > 0.0);
}

// -- Test 14: the subnormal underflow route to a negative duration -----------
TEST_CASE("Trapezoidal: a boundary velocity whose square underflows is rejected",
          "[traj][trapezoidal]")
{
    using traj_t = ctrlpp::trapezoidal_trajectory<double>;

    // Both eq. (3.14)'s feasibility test and the triangular peak
    // sqrt(a h + (v0^2 + v1^2) / 2) are built from the SQUARES of the boundary
    // velocities. Below the square root of the smallest normal value a boundary
    // velocity squares into the subnormal range and starts losing digits; below
    // the square root of the smallest subnormal it squares to exactly zero, the
    // feasibility test reads as satisfied on a command that does not satisfy it,
    // and the peak lands below the boundary velocity it is analytically bounded
    // by. T_a = (v_peak - v0) / a is then NEGATIVE, and evaluate() would reach
    // std::clamp(t, 0, T) with the lower bound above the upper one, which is
    // undefined behavior rather than an odd clamp. Both boundaries are
    // properties of the scalar type, so both are spelled as ones.
    auto const squares_to_subnormal = std::sqrt(std::numeric_limits<double>::min()) / 2.0;
    auto const squares_to_zero = std::sqrt(std::numeric_limits<double>::denorm_min()) / 2.0;
    REQUIRE(squares_to_zero > 0.0);
    REQUIRE(squares_to_zero * squares_to_zero == 0.0);
    REQUIRE(squares_to_subnormal > squares_to_zero);

    for (double const below : {squares_to_subnormal, squares_to_zero}) {
        CAPTURE(below);
        auto const rejected = traj_t::create(
            {.q0 = -1.0, .q1 = -1.0, .v_max = 1e-6, .a_max = 1e-6, .v0 = below, .v1 = 0.0});
        REQUIRE_FALSE(rejected.has_value());
        REQUIRE(rejected.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);
    }

    // The exact artifact the fuzz oracle surfaced once a negative duration
    // became an abort condition.
    auto const artifact = traj_t::create({.q0 = -1.718333114002038e-93,
                                          .q1 = -1.718333114002038e-93,
                                          .v_max = 1e-6,
                                          .a_max = 1e-6,
                                          .v0 = 3.664e-312,
                                          .v1 = 0.0});
    REQUIRE_FALSE(artifact.has_value());
    REQUIRE(artifact.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);

    // The same command with a displacement large enough for the algebra to hold
    // is served, so the rejection is confined to the regime that motivates it.
    auto const served = traj_t::create(
        {.q0 = 0.0, .q1 = 1.0, .v_max = 1e-6, .a_max = 1e-6, .v0 = squares_to_zero, .v1 = 0.0});
    REQUIRE(served.has_value());
    REQUIRE(served.value().duration() > 0.0);
}

// -- Test 15: the input domain of the limits and the boundary values ---------
TEST_CASE("Trapezoidal: out-of-domain limits and non-finite inputs are rejected",
          "[traj][trapezoidal]")
{
    using traj_t = ctrlpp::trapezoidal_trajectory<double>;
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    // The velocity limit is the cruise velocity the profile holds and the
    // divisor of the cruise duration, so its domain is finite and strictly
    // positive. A negative one previously produced a negative acceleration
    // phase, which is the same undefined clamp by another route.
    for (double const v_max : {0.0, -1.0, nan, inf}) {
        auto const result = traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = v_max, .a_max = 1.0});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_velocity_limit);
    }

    // The acceleration divides both ramp durations.
    for (double const a_max : {0.0, -1.0, nan, inf}) {
        auto const result = traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = a_max});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_acceleration_limit);
    }

    for (double const bad : {nan, inf, -inf}) {
        REQUIRE(traj_t::create({.q0 = bad, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0}).error()
                == ctrlpp::trajectory_error::non_finite_input);
        REQUIRE(traj_t::create({.q0 = 0.0, .q1 = bad, .v_max = 1.0, .a_max = 1.0}).error()
                == ctrlpp::trajectory_error::non_finite_input);
        REQUIRE(traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .v0 = bad})
                    .error()
                == ctrlpp::trajectory_error::non_finite_input);
        REQUIRE(traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .v1 = bad})
                    .error()
                == ctrlpp::trajectory_error::non_finite_input);
    }
}

// -- Test 16: a profile that does not fit the scalar type --------------------
TEST_CASE("Trapezoidal: a duration outside the representable range is rejected",
          "[traj][trapezoidal]")
{
    using traj_t = ctrlpp::trapezoidal_trajectory<double>;

    // The cruise duration is the residual displacement over the cruise velocity.
    // A displacement near the top of the range divided by a velocity limit near
    // the bottom of it overflows, and an infinite duration is not a profile: it
    // is a command this scalar type cannot describe.
    auto const overflowed =
        traj_t::create({.q0 = 0.0, .q1 = 1e308, .v_max = 1e-300, .a_max = 1.0});
    REQUIRE_FALSE(overflowed.has_value());
    REQUIRE(overflowed.error() == ctrlpp::trajectory_error::unrepresentable_duration);

    // The same command at a limit that keeps the quotient inside the range is
    // served, so the rejection tracks representability rather than magnitude.
    auto const served = traj_t::create({.q0 = 0.0, .q1 = 1e308, .v_max = 1e6, .a_max = 1.0});
    REQUIRE(served.has_value());
    REQUIRE(std::isfinite(served.value().duration()));
}
