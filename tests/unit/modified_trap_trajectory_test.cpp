/// @brief Tests for modified trapezoidal velocity profile.
///
/// Verifies the six-phase modified trapezoidal profile from B&M Sec 3.7.
/// The profile uses cycloidal acceleration ramps instead of step changes,
/// producing smoother motion with reduced vibration excitation.
///
/// @cite biagiotti2009 -- Sec. 3.7, p.119-122

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "ctrlpp/trajectory/modified_trap_trajectory.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("Modified trapezoidal: endpoint accuracy", "[traj][modified_trap]")
{
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(2.0);

    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
}

TEST_CASE("Modified trapezoidal: symmetry about T/2", "[traj][modified_trap]")
{
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    auto const mid = traj.evaluate(1.0);

    CHECK_THAT(mid.position[0], WithinAbs(5.0, 1e-12));
}

TEST_CASE("Modified trapezoidal: zero endpoint velocities", "[traj][modified_trap]")
{
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(2.0);

    CHECK_THAT(p0.velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Modified trapezoidal: peak velocity == 2*h/T", "[traj][modified_trap]")
{
    // B&M Sec 3.7: max velocity = 2*h/T
    // h=10, T=2 -> max_vel = 2*10/2 = 10.0
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    CHECK_THAT(traj.peak_velocity(), WithinAbs(10.0, 1e-10));

    // Verify at T/2 where velocity peaks
    auto const mid = traj.evaluate(1.0);
    CHECK_THAT(mid.velocity[0], WithinAbs(10.0, 1e-10));
}

TEST_CASE("Modified trapezoidal: peak acceleration == 4.888*h/T^2", "[traj][modified_trap]")
{
    // B&M Sec 3.7: max acc = 2*h*(2+pi)/(pi*T^2) ~= 4.888*h/T^2
    // h=10, T=2 -> max_acc ~= 4.888*10/4 = 12.22
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    constexpr double pi = 3.14159265358979323846;
    // a_max = hp * 8*pi/T^2, where hp = h/(2+pi)
    // = h * 8*pi / ((2+pi) * T^2)
    double const expected_max_acc = 10.0 * 8.0 * pi / ((2.0 + pi) * 4.0);

    // Sample acceleration at T/4 where it should be at maximum (constant accel region)
    auto const pt = traj.evaluate(0.5);
    CHECK_THAT(pt.acceleration[0], WithinAbs(expected_max_acc, 1e-10));
}

TEST_CASE("Modified trapezoidal: phase boundary continuity", "[traj][modified_trap]")
{
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    double const T = 2.0;
    // Phase boundaries at T/8, T/4 (=2T/8), 3T/8, T/2, 5T/8, 3T/4 (=6T/8), 7T/8
    double const boundaries[] = {T / 8.0, T / 4.0, 3.0 * T / 8.0, T / 2.0,
                                  5.0 * T / 8.0, 3.0 * T / 4.0, 7.0 * T / 8.0};
    double const eps = 1e-13;

    for (double b : boundaries) {
        auto const left = traj.evaluate(b - eps);
        auto const right = traj.evaluate(b + eps);

        CAPTURE(b);
        CHECK_THAT(left.position[0], WithinAbs(right.position[0], 1e-10));
        CHECK_THAT(left.velocity[0], WithinAbs(right.velocity[0], 1e-8));
    }
}

TEST_CASE("Modified trapezoidal: negative displacement", "[traj][modified_trap]")
{
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 10.0, .q1 = 0.0, .T = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(2.0);
    auto const mid = traj.evaluate(1.0);

    CHECK_THAT(p0.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(mid.position[0], WithinAbs(5.0, 1e-12));

    // Velocity should be negative
    CHECK(mid.velocity[0] < 0.0);
}

TEST_CASE("Modified trapezoidal: satisfies trajectory_segment concept", "[traj][modified_trap]")
{
    static_assert(ctrlpp::trajectory_segment<
        ctrlpp::modified_trap_trajectory<double>, double, 1>);
    SUCCEED("Concept satisfied");
}

TEST_CASE("Modified trapezoidal: introspection", "[traj][modified_trap]")
{
    ctrlpp::modified_trap_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    CHECK_THAT(traj.duration(), WithinAbs(2.0, 1e-15));
    CHECK(traj.peak_velocity() > 0.0);
}
