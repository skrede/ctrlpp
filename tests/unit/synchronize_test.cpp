#include "ctrlpp/traj/synchronize.h"
#include "ctrlpp/traj/double_s_trajectory.h"
#include "ctrlpp/traj/trapezoidal_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;

// --------------------------------------------------------------------------
// Test 1: Three trapezoidal axes synchronize to slowest
// --------------------------------------------------------------------------
TEST_CASE("synchronize: 3 trapezoidal axes equal duration", "[traj][sync]")
{
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax3({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});

    auto const max_dur = std::max({ax1.duration(), ax2.duration(), ax3.duration()});

    ctrlpp::synchronize(ax1, ax2, ax3);

    REQUIRE_THAT(ax1.duration(), WithinAbs(max_dur, 1e-10));
    REQUIRE_THAT(ax2.duration(), WithinAbs(max_dur, 1e-10));
    REQUIRE_THAT(ax3.duration(), WithinAbs(max_dur, 1e-10));
}

// --------------------------------------------------------------------------
// Test 2: Two double-S axes synchronize to slowest
// --------------------------------------------------------------------------
TEST_CASE("synchronize: 2 double-S axes equal duration", "[traj][sync]")
{
    ctrlpp::double_s_trajectory<double> ax1(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});
    ctrlpp::double_s_trajectory<double> ax2(
        {.q0 = 0.0, .q1 = 3.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});

    auto const max_dur = std::max(ax1.duration(), ax2.duration());

    ctrlpp::synchronize(ax1, ax2);

    REQUIRE_THAT(ax1.duration(), WithinAbs(max_dur, 1e-10));
    REQUIRE_THAT(ax2.duration(), WithinAbs(max_dur, 1e-10));
}

// --------------------------------------------------------------------------
// Test 3: Heterogeneous mix -- trapezoidal + double-S
// --------------------------------------------------------------------------
TEST_CASE("synchronize: heterogeneous trapezoidal + double-S", "[traj][sync]")
{
    ctrlpp::trapezoidal_trajectory<double> trap({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::double_s_trajectory<double> ds(
        {.q0 = 0.0, .q1 = 3.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});

    auto const max_dur = std::max(trap.duration(), ds.duration());

    ctrlpp::synchronize(trap, ds);

    REQUIRE_THAT(trap.duration(), WithinAbs(max_dur, 1e-10));
    REQUIRE_THAT(ds.duration(), WithinAbs(max_dur, 1e-10));
}

// --------------------------------------------------------------------------
// Test 4: Post-sync endpoints are exact
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync endpoints preserved", "[traj][sync]")
{
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 1.0, .q1 = 11.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 2.0, .q1 = 22.0, .v_max = 5.0, .a_max = 10.0});

    ctrlpp::synchronize(ax1, ax2);

    auto const p1_start = ax1.evaluate(0.0);
    auto const p1_end = ax1.evaluate(ax1.duration());
    REQUIRE_THAT(p1_start.position[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(p1_end.position[0], WithinAbs(11.0, 1e-10));

    auto const p2_start = ax2.evaluate(0.0);
    auto const p2_end = ax2.evaluate(ax2.duration());
    REQUIRE_THAT(p2_start.position[0], WithinAbs(2.0, 1e-10));
    REQUIRE_THAT(p2_end.position[0], WithinAbs(22.0, 1e-10));
}

// --------------------------------------------------------------------------
// Test 5: Post-sync velocity never exceeds v_max
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync velocity within v_max", "[traj][sync]")
{
    double constexpr v_max = 5.0;
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 0.0, .q1 = 10.0, .v_max = v_max, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 0.0, .q1 = 20.0, .v_max = v_max, .a_max = 10.0});

    ctrlpp::synchronize(ax1, ax2);

    double constexpr eps = 1e-10;
    for (int i = 0; i <= 1000; ++i) {
        double const t = ax1.duration() * static_cast<double>(i) / 1000.0;
        auto const pt = ax1.evaluate(t);
        REQUIRE(std::abs(pt.velocity[0]) <= v_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 6: Post-sync acceleration never exceeds a_max
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync acceleration within a_max", "[traj][sync]")
{
    double constexpr a_max = 10.0;
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = a_max});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = a_max});

    ctrlpp::synchronize(ax1, ax2);

    double constexpr eps = 1e-10;
    for (int i = 0; i <= 1000; ++i) {
        double const t = ax1.duration() * static_cast<double>(i) / 1000.0;
        auto const pt = ax1.evaluate(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= a_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 7: Single axis is a no-op
// --------------------------------------------------------------------------
TEST_CASE("synchronize: single axis no-op", "[traj][sync]")
{
    ctrlpp::trapezoidal_trajectory<double> ax({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    auto const dur_before = ax.duration();

    ctrlpp::synchronize(ax);

    REQUIRE_THAT(ax.duration(), WithinAbs(dur_before, 1e-14));
}

// --------------------------------------------------------------------------
// Test 8: Vector overload with trapezoidal axes
// --------------------------------------------------------------------------
TEST_CASE("synchronize: vector overload", "[traj][sync]")
{
    std::vector<ctrlpp::trapezoidal_trajectory<double>> axes;
    axes.emplace_back(ctrlpp::trapezoidal_trajectory<double>::config{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    axes.emplace_back(ctrlpp::trapezoidal_trajectory<double>::config{
        .q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});
    axes.emplace_back(ctrlpp::trapezoidal_trajectory<double>::config{
        .q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0});

    double max_dur = 0.0;
    for (auto const& ax : axes) {
        max_dur = std::max(max_dur, ax.duration());
    }

    ctrlpp::synchronize(axes);

    for (auto const& ax : axes) {
        REQUIRE_THAT(ax.duration(), WithinAbs(max_dur, 1e-10));
    }
}
