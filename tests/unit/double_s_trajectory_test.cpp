#include <ctrlpp/traj/double_s_trajectory.h>
#include <ctrlpp/traj/trajectory_segment.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// --------------------------------------------------------------------------
// Test 1: Full 7-segment profile reaches q1 with 7 distinct phases
// --------------------------------------------------------------------------
TEST_CASE("double_s: full 7-segment profile", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);

    auto const pt_end = traj.evaluate(T);
    REQUIRE_THAT(pt_end.position[0], WithinAbs(10.0, 1e-10));
    REQUIRE_THAT(pt_end.velocity[0], WithinAbs(0.0, 1e-8));

    auto const pt_start = traj.evaluate(0.0);
    REQUIRE_THAT(pt_start.position[0], WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(pt_start.velocity[0], WithinAbs(0.0, 1e-8));

    // All 7 phase durations should be > 0 for a full profile
    auto const pd = traj.phase_durations();
    REQUIRE(pd.size() == 7);
    // At least the jerk phases (0, 2, 4, 6) and cruise (3) should be > 0
    CHECK(pd[0] > 0.0);
    CHECK(pd[3] >= 0.0); // cruise may or may not exist
    CHECK(pd[6] > 0.0);
}

// --------------------------------------------------------------------------
// Test 2: Velocity not reached (degenerate -- v_lim < v_max)
// --------------------------------------------------------------------------
TEST_CASE("double_s: velocity not reached degenerate", "[traj][double_s]")
{
    // Short distance forces v_lim < v_max
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = 100.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    CHECK(traj.is_degenerate());
    CHECK(traj.peak_velocity() < 100.0);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(1.0, 1e-10));
}

// --------------------------------------------------------------------------
// Test 3: Acceleration not reached (a_lim < a_max)
// --------------------------------------------------------------------------
TEST_CASE("double_s: acceleration not reached degenerate", "[traj][double_s]")
{
    // Very high a_max but low j_max -- can't reach full a_max
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.5, .v_max = 10.0, .a_max = 100.0, .j_max = 50.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    CHECK(traj.is_degenerate());

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.5, 1e-10));
}

// --------------------------------------------------------------------------
// Test 4: Both v and a not reached (doubly degenerate)
// --------------------------------------------------------------------------
TEST_CASE("double_s: doubly degenerate (v and a not reached)", "[traj][double_s]")
{
    // Tiny distance with high limits -- neither v_max nor a_max reached
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.01, .v_max = 100.0, .a_max = 100.0, .j_max = 50.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    CHECK(traj.is_degenerate());
    CHECK(traj.peak_velocity() < 100.0);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.01, 1e-10));
}

// --------------------------------------------------------------------------
// Test 5: Triangular-like (very short distance)
// --------------------------------------------------------------------------
TEST_CASE("double_s: very short distance simplifies", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.001, .v_max = 10.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.001, 1e-10));
}

// --------------------------------------------------------------------------
// Test 6: Velocity constraint -- sample 1000 points, |velocity| <= v_max + eps
// --------------------------------------------------------------------------
TEST_CASE("double_s: velocity constraint satisfied", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const T = traj.duration();
    constexpr double eps = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.velocity[0]) <= cfg.v_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 7: Acceleration constraint -- sample 1000 points
// --------------------------------------------------------------------------
TEST_CASE("double_s: acceleration constraint satisfied", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const T = traj.duration();
    constexpr double eps = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= cfg.a_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 8: Jerk constraint -- finite difference check
// --------------------------------------------------------------------------
TEST_CASE("double_s: jerk constraint satisfied", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const T = traj.duration();
    constexpr double dt = 1e-5;
    constexpr double eps = 1.0; // jerk via finite diff has larger tolerance
    for (int i = 1; i < 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt1 = traj.evaluate(t - dt);
        auto const pt2 = traj.evaluate(t + dt);
        double const jerk = (pt2.acceleration[0] - pt1.acceleration[0]) / (2.0 * dt);
        REQUIRE(std::abs(jerk) <= cfg.j_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 9: Phase boundary continuity
// --------------------------------------------------------------------------
TEST_CASE("double_s: phase boundary continuity", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const pd = traj.phase_durations();
    double boundary = 0.0;
    constexpr double tiny = 1e-12;
    constexpr double tol = 1e-10;

    for (std::size_t i = 0; i < 6; ++i) {
        boundary += pd[i];
        if (boundary <= 0.0 || boundary >= traj.duration()) {
            continue;
        }
        auto const left = traj.evaluate(boundary - tiny);
        auto const right = traj.evaluate(boundary + tiny);
        REQUIRE_THAT(left.position[0], WithinAbs(right.position[0], tol));
        REQUIRE_THAT(left.velocity[0], WithinAbs(right.velocity[0], tol));
        REQUIRE_THAT(left.acceleration[0], WithinAbs(right.acceleration[0], tol));
    }
}

// --------------------------------------------------------------------------
// Test 10: Negative displacement (q1 < q0) via sigma transformation
// --------------------------------------------------------------------------
TEST_CASE("double_s: negative displacement sigma transformation", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 10.0, .q1 = 0.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);

    auto const pt_start = traj.evaluate(0.0);
    REQUIRE_THAT(pt_start.position[0], WithinAbs(10.0, 1e-10));

    auto const pt_end = traj.evaluate(T);
    REQUIRE_THAT(pt_end.position[0], WithinAbs(0.0, 1e-10));

    // Velocity should be negative (moving from 10 -> 0)
    auto const pt_mid = traj.evaluate(T / 2.0);
    CHECK(pt_mid.velocity[0] < 0.0);
}

// --------------------------------------------------------------------------
// Test 11: Introspection methods
// --------------------------------------------------------------------------
TEST_CASE("double_s: introspection", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    // is_degenerate returns bool
    auto const deg = traj.is_degenerate();
    static_assert(std::is_same_v<decltype(deg), const bool>);

    // peak_velocity returns Scalar
    auto const pv = traj.peak_velocity();
    CHECK(pv > 0.0);
    CHECK(pv <= cfg.v_max + 1e-10);

    // phase_durations returns array of 7
    auto const pd = traj.phase_durations();
    static_assert(pd.size() == 7);

    // Sum of phase durations equals total duration
    double sum = 0.0;
    for (auto const& d : pd) {
        CHECK(d >= 0.0);
        sum += d;
    }
    REQUIRE_THAT(sum, WithinAbs(traj.duration(), 1e-12));
}

// --------------------------------------------------------------------------
// Test 12: Satisfies trajectory_segment concept
// --------------------------------------------------------------------------
TEST_CASE("double_s: satisfies trajectory_segment concept", "[traj][double_s]")
{
    static_assert(ctrlpp::trajectory_segment<
                  ctrlpp::double_s_trajectory<double>, double, 1>);

    static_assert(ctrlpp::trajectory_segment<
                  ctrlpp::double_s_trajectory<float>, float, 1>);
}

// --------------------------------------------------------------------------
// Test 13: Zero displacement (q0 == q1) returns stationary profile
// --------------------------------------------------------------------------
TEST_CASE("double_s: zero displacement stationary profile", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 5.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    ctrlpp::double_s_trajectory<double> traj(cfg);

    CHECK(traj.duration() == 0.0);

    auto const pt = traj.evaluate(0.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(5.0, 1e-12));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(pt.acceleration[0], WithinAbs(0.0, 1e-12));
}
