/// @brief Tests for modified sinusoidal velocity profile.
///
/// Verifies the three-region modified sinusoidal profile from B&M Sec 3.8.
/// The profile uses a sinusoidal middle region with cycloidal start/end
/// transitions, producing very smooth acceleration with minimal harmonics.
///
/// @cite biagiotti2009 -- Sec. 3.8, p.124-127

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "ctrlpp/trajectory/modified_sin_trajectory.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("Modified sinusoidal: endpoint accuracy", "[traj][modified_sin]")
{
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(2.0);

    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
}

TEST_CASE("Modified sinusoidal: symmetry about T/2", "[traj][modified_sin]")
{
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    auto const mid = traj.evaluate(1.0);

    CHECK_THAT(mid.position[0], WithinAbs(5.0, 1e-12));
}

TEST_CASE("Modified sinusoidal: zero endpoint velocities", "[traj][modified_sin]")
{
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(2.0);

    CHECK_THAT(p0.velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Modified sinusoidal: peak velocity ~= 1.76*h/T", "[traj][modified_sin]")
{
    // B&M Sec 3.8: max velocity = 4*pi*h / ((pi+4)*T) ~= 1.7596*h/T
    // h=10, T=2 -> max_vel ~= 8.798
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    constexpr double pi = 3.14159265358979323846;
    double const expected_peak_vel = 4.0 * pi * 10.0 / ((pi + 4.0) * 2.0);

    CHECK_THAT(traj.peak_velocity(), WithinAbs(expected_peak_vel, 1e-10));

    // At T/2, velocity should be at its peak
    auto const mid = traj.evaluate(1.0);
    CHECK_THAT(mid.velocity[0], WithinAbs(expected_peak_vel, 1e-10));
}

TEST_CASE("Modified sinusoidal: peak acceleration ~= 5.528*h/T^2", "[traj][modified_sin]")
{
    // B&M Sec 3.8: max acc = 4*pi^2*h / ((pi+4)*T^2) ~= 5.528*h/T^2
    // For h=10, T=2: max_acc ~= 5.528*10/4 = 13.82
    // Peak acceleration occurs at T/8 (boundary between cycloidal and sinusoidal regions)
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    constexpr double pi = 3.14159265358979323846;
    double const expected_max_acc = 4.0 * pi * pi * 10.0 / ((pi + 4.0) * 4.0);

    // Sample at T/8 where acceleration is at its peak
    auto const pt = traj.evaluate(0.25);
    CHECK_THAT(pt.acceleration[0], WithinAbs(expected_max_acc, 1e-8));
}

TEST_CASE("Modified sinusoidal: phase boundary continuity", "[traj][modified_sin]")
{
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    double const T = 2.0;
    // Phase boundaries at T/8, 7T/8, T/2
    double const boundaries[] = {T / 8.0, T / 2.0, 7.0 * T / 8.0};
    double const eps = 1e-13;

    for (double b : boundaries) {
        auto const left = traj.evaluate(b - eps);
        auto const right = traj.evaluate(b + eps);

        CAPTURE(b);
        CHECK_THAT(left.position[0], WithinAbs(right.position[0], 1e-10));
        CHECK_THAT(left.velocity[0], WithinAbs(right.velocity[0], 1e-8));
    }
}

TEST_CASE("Modified sinusoidal: negative displacement", "[traj][modified_sin]")
{
    ctrlpp::modified_sin_trajectory<double> traj({.q0 = 10.0, .q1 = 0.0, .T = 2.0});

    auto const p0 = traj.evaluate(0.0);
    auto const pT = traj.evaluate(2.0);
    auto const mid = traj.evaluate(1.0);

    CHECK_THAT(p0.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(mid.position[0], WithinAbs(5.0, 1e-12));

    // Velocity should be negative
    CHECK(mid.velocity[0] < 0.0);
}

TEST_CASE("Modified sinusoidal: satisfies trajectory_segment concept", "[traj][modified_sin]")
{
    static_assert(ctrlpp::trajectory_segment<
        ctrlpp::modified_sin_trajectory<double>, double, 1>);
    SUCCEED("Concept satisfied");
}
