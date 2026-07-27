#include "ctrlpp/trajectory/synchronize.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <algorithm>

using Catch::Matchers::WithinAbs;

namespace
{

// Displacement swept by a profile, integrated from the velocity it reports over
// its own phase segments. A position sample at the end of the move is not usable
// as an oracle here: both profiles write their final segment as an offset
// backwards from the commanded displacement, so the end position comes back
// correct even when the traversal is not. Panels are laid inside each segment
// because a panel straddling a phase kink carries a truncation error that would
// manufacture violations on correct profiles.
template <typename Profile>
auto swept_displacement(Profile const& profile) -> double
{
    constexpr int panels_per_segment = 8;

    double integral = 0.0;
    double t_start = 0.0;
    for (auto const& segment : profile.phase_durations()) {
        double const length = static_cast<double>(segment);
        if (!(length > 0.0)) {
            continue;
        }
        double const dt = length / static_cast<double>(panels_per_segment);
        for (int p = 0; p < panels_per_segment; ++p) {
            double const a = t_start + static_cast<double>(p) * dt;
            double const b = a + dt;
            double const m = a + 0.5 * dt;
            integral += (dt / 6.0)
                * (profile.evaluate(a).velocity(0) + 4.0 * profile.evaluate(m).velocity(0)
                   + profile.evaluate(b).velocity(0));
        }
        t_start += length;
    }
    return integral;
}

// Velocity a step inside the end of the move, where the step is the profile's
// own final nonempty segment. Nothing anchors the velocity the way the final
// position segment is anchored to the commanded displacement, so this is the
// half of the traversal oracle that a position sample cannot supply.
template <typename Profile>
auto terminal_velocity(Profile const& profile) -> double
{
    double total = 0.0;
    double last = 0.0;
    for (auto const& segment : profile.phase_durations()) {
        double const length = static_cast<double>(segment);
        if (length > 0.0) {
            total += length;
            last = length;
        }
    }
    if (!(last > 0.0)) {
        return 0.0;
    }
    return static_cast<double>(profile.evaluate(total - last).velocity(0));
}

/// Bound on how far the velocity a step inside the end may sit from the
/// commanded final velocity: the acceleration limit times that step, plus one
/// sample's worth of rounding at the velocity scale and the backward-time term
/// the deceleration branch carries (it recovers its local time by subtracting
/// from the duration, which is representable only to one unit in the last place
/// of it, and the velocity slews at up to the acceleration limit there).
auto terminal_velocity_bound(double a_max, double v_scale, double duration, double step) -> double
{
    constexpr int rounding_ops_per_sample = 8;
    return a_max * step
           + static_cast<double>(rounding_ops_per_sample) * std::numeric_limits<double>::epsilon()
                 * (v_scale + a_max * duration);
}

template <typename Profile>
auto final_segment(Profile const& profile) -> double
{
    double last = 0.0;
    for (auto const& segment : profile.phase_durations()) {
        double const length = static_cast<double>(segment);
        if (length > 0.0) {
            last = length;
        }
    }
    return last;
}

/// A profile's duration together with a trace sampled at fixed interior
/// fractions of it, captured so that a later capture can be compared to it bit
/// for bit. Exact comparison is the point: a rejected synchronization must not
/// have moved any axis at all, and "close enough" would not say that.
struct axis_snapshot
{
    static constexpr int samples = 9;

    double duration{};
    std::array<double, samples> position{};
    std::array<double, samples> velocity{};
    std::array<double, samples> acceleration{};
};

template <typename Profile>
auto snapshot(Profile const& profile) -> axis_snapshot
{
    axis_snapshot out{};
    out.duration = static_cast<double>(profile.duration());
    for (int i = 0; i < axis_snapshot::samples; ++i) {
        double const fraction = static_cast<double>(i + 1) / static_cast<double>(axis_snapshot::samples + 1);
        auto const point = profile.evaluate(out.duration * fraction);
        out.position[static_cast<std::size_t>(i)] = static_cast<double>(point.position(0));
        out.velocity[static_cast<std::size_t>(i)] = static_cast<double>(point.velocity(0));
        out.acceleration[static_cast<std::size_t>(i)] = static_cast<double>(point.acceleration(0));
    }
    return out;
}

auto identical(axis_snapshot const& lhs, axis_snapshot const& rhs) -> bool
{
    if (lhs.duration != rhs.duration) {
        return false;
    }
    for (int i = 0; i < axis_snapshot::samples; ++i) {
        auto const k = static_cast<std::size_t>(i);
        if (lhs.position[k] != rhs.position[k] || lhs.velocity[k] != rhs.velocity[k]
            || lhs.acceleration[k] != rhs.acceleration[k]) {
            return false;
        }
    }
    return true;
}

// An axis whose commanded displacement sits below the distance its two boundary
// velocities already sweep between them. The reachable durations of such a
// profile are bounded above, so it refuses any target beyond that supremum. Its
// own duration is about 0.47 s and its supremum about 0.54 s.
auto bounded_axis() -> ctrlpp::trapezoidal_trajectory<double>
{
    return ctrlpp::trapezoidal_trajectory<double>(
        {.q0 = 0.0, .q1 = 1.0, .v_max = 3.0, .a_max = 1.0, .v0 = 2.0, .v1 = 2.0});
}

// A slower axis with room to spare, so it sets the synchronization target.
auto slow_axis() -> ctrlpp::trapezoidal_trajectory<double>
{
    return ctrlpp::trapezoidal_trajectory<double>(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
}

}

// --------------------------------------------------------------------------
// Test 1: Three trapezoidal axes synchronize to slowest
// --------------------------------------------------------------------------
TEST_CASE("synchronize: 3 trapezoidal axes equal duration", "[traj][sync]")
{
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax3({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});

    auto const max_dur = std::max({ax1.duration(), ax2.duration(), ax3.duration()});

    REQUIRE(ctrlpp::synchronize(ax1, ax2, ax3).has_value());

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

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

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

    REQUIRE(ctrlpp::synchronize(trap, ds).has_value());

    REQUIRE_THAT(trap.duration(), WithinAbs(max_dur, 1e-10));
    REQUIRE_THAT(ds.duration(), WithinAbs(max_dur, 1e-10));
}

// --------------------------------------------------------------------------
// Test 4: Post-sync traversal is preserved
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync traversal preserved", "[traj][sync]")
{
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 1.0, .q1 = 11.0, .v_max = 5.0, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 2.0, .q1 = 22.0, .v_max = 5.0, .a_max = 10.0});

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

    // The start position is a real value the profile has to produce.
    auto const p1_start = ax1.evaluate(0.0);
    auto const p2_start = ax2.evaluate(0.0);
    REQUIRE_THAT(p1_start.position[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(p2_start.position[0], WithinAbs(2.0, 1e-10));

    // The traversal is asserted by integrating the reported velocity, not by
    // sampling the end position, which the final segment reproduces by
    // construction whatever the profile actually does.
    REQUIRE_THAT(swept_displacement(ax1), WithinAbs(10.0, 1e-10));
    REQUIRE_THAT(swept_displacement(ax2), WithinAbs(20.0, 1e-10));

    // Second half of the oracle: both axes were commanded to come to rest, and
    // the velocity a step inside the end is not anchored to anything, so it can
    // report a profile that ends still moving.
    double constexpr a_max = 10.0;
    double constexpr v_max = 5.0;
    for (auto const* axis : {&ax1, &ax2}) {
        double const step = final_segment(*axis);
        REQUIRE(step > 0.0);
        REQUIRE(std::abs(terminal_velocity(*axis))
                <= terminal_velocity_bound(a_max, v_max, axis->duration(), step));
    }
}

// --------------------------------------------------------------------------
// Test 5: Post-sync velocity never exceeds v_max
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync velocity within v_max", "[traj][sync]")
{
    double constexpr v_max = 5.0;
    ctrlpp::trapezoidal_trajectory<double> ax1({.q0 = 0.0, .q1 = 10.0, .v_max = v_max, .a_max = 10.0});
    ctrlpp::trapezoidal_trajectory<double> ax2({.q0 = 0.0, .q1 = 20.0, .v_max = v_max, .a_max = 10.0});

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

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

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

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

    REQUIRE(ctrlpp::synchronize(ax).has_value());

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

    REQUIRE(ctrlpp::synchronize(std::span{axes}).has_value());

    for (auto const& ax : axes) {
        REQUIRE_THAT(ax.duration(), WithinAbs(max_dur, 1e-10));
    }
}

// --------------------------------------------------------------------------
// Test 9: A rejected synchronization leaves every axis bitwise untouched
//         (variadic overload)
// --------------------------------------------------------------------------
TEST_CASE("synchronize: rejection mutates no axis, variadic overload", "[traj][sync]")
{
    auto ax1 = slow_axis();
    auto ax2 = bounded_axis();
    auto ax3 = slow_axis();

    // The target is the slowest axis's duration, which lies past what the
    // bounded axis can ever reach.
    REQUIRE(ax1.duration() > ax2.duration());

    auto const before1 = snapshot(ax1);
    auto const before2 = snapshot(ax2);
    auto const before3 = snapshot(ax3);

    auto const result = ctrlpp::synchronize(ax1, ax2, ax3);
    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::trajectory_error::unreachable_duration);

    REQUIRE(identical(before1, snapshot(ax1)));
    REQUIRE(identical(before2, snapshot(ax2)));
    REQUIRE(identical(before3, snapshot(ax3)));
}

// --------------------------------------------------------------------------
// Test 10: The same, through the contiguous-view overload
// --------------------------------------------------------------------------
TEST_CASE("synchronize: rejection mutates no axis, contiguous-view overload", "[traj][sync]")
{
    std::array<ctrlpp::trapezoidal_trajectory<double>, 3> axes{slow_axis(), bounded_axis(),
                                                               slow_axis()};

    std::array<axis_snapshot, 3> before{};
    for (std::size_t i = 0; i < axes.size(); ++i) {
        before[i] = snapshot(axes[i]);
    }

    auto const result = ctrlpp::synchronize(std::span<ctrlpp::trapezoidal_trajectory<double>>{axes});
    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::trajectory_error::unreachable_duration);

    for (std::size_t i = 0; i < axes.size(); ++i) {
        REQUIRE(identical(before[i], snapshot(axes[i])));
    }
}

// --------------------------------------------------------------------------
// Test 11: The slowest axis takes the success no-op path
// --------------------------------------------------------------------------
TEST_CASE("synchronize: the slowest axis is a success no-op", "[traj][sync]")
{
    auto slowest = slow_axis();
    ctrlpp::trapezoidal_trajectory<double> faster(
        {.q0 = 0.0, .q1 = 2.0, .v_max = 5.0, .a_max = 10.0});
    REQUIRE(slowest.duration() > faster.duration());

    auto const slowest_before = snapshot(slowest);

    REQUIRE(ctrlpp::synchronize(slowest, faster).has_value());

    // Exact equality is safe here, and this is not float-equality fragility: the
    // target duration is a bit-exact copy of this axis's own stored duration, so
    // the retiming compares equal to it and returns without touching anything,
    // rather than taking the shortening rejection or rebuilding the profile.
    REQUIRE(identical(slowest_before, snapshot(slowest)));

    // The other axis was genuinely rebuilt, so its duration is not bit-exact and
    // must not be asserted as such. Twenty-seven chained roundings stand behind
    // it: twelve form the cruise velocity (the linear coefficient, the constant,
    // the discriminant, its square root, and the cancellation-free root
    // selection) and fifteen form the duration recomputed from it. Both boundary
    // velocities are zero here, so the shape residual is a sum of nonnegative
    // terms and carries no cancellation to amplify them.
    constexpr int duration_rounding_ops = 27;
    REQUIRE(std::abs(faster.duration() - slowest.duration())
            <= static_cast<double>(duration_rounding_ops) * std::numeric_limits<double>::epsilon()
                   * slowest.duration());
}

// --------------------------------------------------------------------------
// Test 12: An empty view is a successful no-op
// --------------------------------------------------------------------------
TEST_CASE("synchronize: empty view succeeds", "[traj][sync]")
{
    std::array<ctrlpp::trapezoidal_trajectory<double>, 0> none{};
    REQUIRE(ctrlpp::synchronize(std::span<ctrlpp::trapezoidal_trajectory<double>>{none}).has_value());
}
