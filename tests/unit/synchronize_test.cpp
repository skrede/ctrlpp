// What the oracles in this file decide.
//
// This file is NOT a duration-equality check with nothing behind it. A prior
// phase gave it nine oracles that decide real properties, and they are the
// templates the rest of the file now follows:
//
//  * The traversal is asserted by a kink-aligned Simpson quadrature of the
//    reported velocity, because an end-position sample is reproduced by
//    construction whatever the profile does (the helper below says so).
//  * The terminal velocity a step inside the end is asserted against a DERIVED
//    bound: the acceleration limit times that step plus a counted eight-operation
//    rounding term.
//  * A rejected synchronization is asserted to have mutated nothing, BITWISE,
//    across the duration and nine interior samples of position, velocity and
//    acceleration, through both overloads, and to carry the specific enumerator.
//  * The slowest axis's success no-op is asserted bitwise, and the axis that was
//    genuinely rebuilt against a twenty-seven-operation counted budget with all
//    twenty-seven enumerated in prose.
//
// What this plan changed is the fifteen thresholds that were not derived from
// anything. Each is now either a counted-operation budget at the scale of the
// operand that entered, or an exact comparison where the quantity is exact:
//
//  * The equal-duration sites take a budget counted for the profile family that
//    was actually rebuilt. The two families have DIFFERENT counts and do not
//    share one because it happens to be larger.
//  * Evaluating at time zero returns the commanded start position with no
//    arithmetic intervening, so those sites are exact equality.
//  * The swept-displacement sites reuse the published per-panel budget from the
//    rescale anchor file, for the same quadrature.
//  * The two envelope sites need no tolerance at all. The reported acceleration
//    is a stored magnitude negated, never computed, and the sampled peak velocity
//    is the profile's own reported peak, so both are exact equalities against the
//    limit and against the profile's own report.
//  * The single-axis no-op reuses the file's own bitwise snapshot comparison,
//    because a tolerance there is not merely unprincipled, it is the wrong kind
//    of assertion for an operation that touches nothing.

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
constexpr int panels_per_segment = 8;

/// Per-panel rounding budget for the kink-aligned quadrature below, published at
/// `trajectory_rescale_anchor_test.cpp:81` for the same quadrature and reused
/// rather than recounted. Two fresh velocity samples per panel at up to eight
/// chained rounding operations each, three multiplies and three adds over those
/// samples, and one addition into each of the two running sums.
constexpr int quadrature_rounding_ops_per_panel = 2 * 8 + 6 + 2;

/// Number of panels the quadrature lays for a given profile. The panels sit
/// inside each nonempty segment, never across a kink.
template <typename Profile>
auto quadrature_panels(Profile const& profile) -> int
{
    int segments = 0;
    for (auto const& segment : profile.phase_durations()) {
        if (static_cast<double>(segment) > 0.0) {
            ++segments;
        }
    }
    return segments * panels_per_segment;
}

/// Chained rounding operations behind a trapezoidal axis's duration after a
/// rebuild at a lower cruise velocity. This is the count enumerated in prose at
/// the slowest-axis case at the end of this file -- twelve forming the cruise
/// velocity and fifteen forming the duration recomputed from it -- named here so
/// the equal-duration cases can use the same number without restating the
/// derivation. Realized: below one unit in the last place of the target.
constexpr int trapezoidal_duration_rounding_ops = 27;

/// The same for a double-S axis, counted separately because its closed form is a
/// different one and sharing the trapezoidal number merely because it is larger
/// would state a budget nobody derived.
///
/// With both boundary velocities at rest the retiming solves the scale in one
/// step and rebuilds through `create`, so the chain runs: the scale itself (one);
/// the scaled velocity limit (two), acceleration limit (three) and jerk limit
/// (four); then, inside the rest-to-rest solve, the reciprocal of the scaled
/// acceleration (four) and the jerk ratio (five) forming the quadratic in the
/// reached velocity, its square of the linear coefficient (six), the product of
/// the constant term (five) and their difference (seven), the square root
/// (eight), the sum with the negated linear coefficient (nine) and the division
/// by twice the leading coefficient (ten); then the ramp duration as that
/// velocity over the scaled acceleration (eleven) added to the jerk ratio
/// (twelve); and finally the total as the two symmetric ramps summed (thirteen).
/// The cruise-present branch of the same solve is shorter, at nine, so thirteen
/// is the longer of the two closed forms and bounds both. Every operation is
/// counted whether or not it rounds. Realized: below one unit in the last place
/// of the target.
constexpr int double_s_duration_rounding_ops = 13;

template <typename Profile>
auto swept_displacement(Profile const& profile) -> double
{
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

// create() is the only construction path on either profile type and it is
// fallible, so every axis a test uses is built through one of these helpers,
// which assert the command was realizable rather than letting a rejection pass
// as an axis that never moves.
auto trapezoidal_axis(ctrlpp::trapezoidal_trajectory<double>::config const& cfg)
    -> ctrlpp::trapezoidal_trajectory<double>
{
    auto created = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

auto double_s_axis(ctrlpp::double_s_trajectory<double>::config const& cfg)
    -> ctrlpp::double_s_trajectory<double>
{
    auto created = ctrlpp::double_s_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

// An axis whose commanded displacement sits below the distance its two boundary
// velocities already sweep between them. The reachable durations of such a
// profile are bounded above, so it refuses any target beyond that supremum. Its
// own duration is about 0.47 s and its supremum about 0.54 s.
auto bounded_axis() -> ctrlpp::trapezoidal_trajectory<double>
{
    return trapezoidal_axis({.q0 = 0.0, .q1 = 1.0, .v_max = 3.0, .a_max = 1.0, .v0 = 2.0, .v1 = 2.0});
}

// A slower axis with room to spare, so it sets the synchronization target.
auto slow_axis() -> ctrlpp::trapezoidal_trajectory<double>
{
    return trapezoidal_axis({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
}

}

// --------------------------------------------------------------------------
// Test 1: Three trapezoidal axes synchronize to slowest
// --------------------------------------------------------------------------
TEST_CASE("synchronize: 3 trapezoidal axes equal duration", "[traj][sync]")
{
    auto ax1 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    auto ax2 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0});
    auto ax3 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});

    auto const max_dur = std::max({ax1.duration(), ax2.duration(), ax3.duration()});

    REQUIRE(ctrlpp::synchronize(ax1, ax2, ax3).has_value());

    // Each retimed axis must reach the target to the precision its own rebuild
    // carries, which the file already counts: twenty-seven operations for a
    // trapezoidal rebuild from rest, times the scalar epsilon, times the target
    // as the scale. The absolute constant this replaces was about four hundred
    // and fifty thousand of those epsilons at this scale, four decades looser
    // than the arithmetic permits.
    double const budget = trapezoidal_duration_rounding_ops
                          * std::numeric_limits<double>::epsilon() * max_dur;
    REQUIRE_THAT(ax1.duration(), WithinAbs(max_dur, budget));
    REQUIRE_THAT(ax2.duration(), WithinAbs(max_dur, budget));
    REQUIRE_THAT(ax3.duration(), WithinAbs(max_dur, budget));
}

// --------------------------------------------------------------------------
// Test 2: Two double-S axes synchronize to slowest
// --------------------------------------------------------------------------
TEST_CASE("synchronize: 2 double-S axes equal duration", "[traj][sync]")
{
    auto ax1 = double_s_axis(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});
    auto ax2 = double_s_axis(
        {.q0 = 0.0, .q1 = 3.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});

    auto const max_dur = std::max(ax1.duration(), ax2.duration());

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

    // The double-S rebuild is a different closed form from the trapezoidal one,
    // so it carries its own count. Thirteen, derived above.
    double const budget = double_s_duration_rounding_ops
                          * std::numeric_limits<double>::epsilon() * max_dur;
    REQUIRE_THAT(ax1.duration(), WithinAbs(max_dur, budget));
    REQUIRE_THAT(ax2.duration(), WithinAbs(max_dur, budget));
}

// --------------------------------------------------------------------------
// Test 3: Heterogeneous mix -- trapezoidal + double-S
// --------------------------------------------------------------------------
TEST_CASE("synchronize: heterogeneous trapezoidal + double-S", "[traj][sync]")
{
    auto trap = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});
    auto ds = double_s_axis(
        {.q0 = 0.0, .q1 = 3.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});

    auto const max_dur = std::max(trap.duration(), ds.duration());

    REQUIRE(ctrlpp::synchronize(trap, ds).has_value());

    // The mixed case is where sharing one number would be wrong: each axis is
    // held to the budget of the family that rebuilt it, not to the larger of the
    // two.
    double const eps = std::numeric_limits<double>::epsilon();
    REQUIRE_THAT(trap.duration(),
                 WithinAbs(max_dur, trapezoidal_duration_rounding_ops * eps * max_dur));
    REQUIRE_THAT(ds.duration(),
                 WithinAbs(max_dur, double_s_duration_rounding_ops * eps * max_dur));
}

// --------------------------------------------------------------------------
// Test 4: Post-sync traversal is preserved
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync traversal preserved", "[traj][sync]")
{
    auto ax1 = trapezoidal_axis(
        {.q0 = 1.0, .q1 = 11.0, .v_max = 5.0, .a_max = 10.0});
    auto ax2 = trapezoidal_axis(
        {.q0 = 2.0, .q1 = 22.0, .v_max = 5.0, .a_max = 10.0});

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

    // The start position is a real value the profile has to produce, and it
    // produces it EXACTLY. Confirmed from the evaluation path rather than
    // assumed: at time zero the acceleration branch forms the position as the
    // initial velocity times zero plus half the acceleration times zero squared,
    // every product an exact product of zero, and the sign transform then adds
    // that exact zero to the commanded start. No arithmetic intervenes, so a
    // tolerance would admit a start the profile cannot reach.
    auto const p1_start = ax1.evaluate(0.0);
    auto const p2_start = ax2.evaluate(0.0);
    REQUIRE(p1_start.position[0] == 1.0);
    REQUIRE(p2_start.position[0] == 2.0);

    // The traversal is asserted by integrating the reported velocity, not by
    // sampling the end position, which the final segment reproduces by
    // construction whatever the profile actually does.
    //
    // Its budget is the per-panel count published for this same quadrature in the
    // rescale anchor file, times the panels this profile carries, times the
    // scalar epsilon, times the commanded displacement -- the operand whose
    // accumulated absolute area the roundings sit on.
    double const eps = std::numeric_limits<double>::epsilon();
    double const budget1 = quadrature_rounding_ops_per_panel * quadrature_panels(ax1) * eps * 10.0;
    double const budget2 = quadrature_rounding_ops_per_panel * quadrature_panels(ax2) * eps * 20.0;
    REQUIRE_THAT(swept_displacement(ax1), WithinAbs(10.0, budget1));
    REQUIRE_THAT(swept_displacement(ax2), WithinAbs(20.0, budget2));

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
    auto ax1 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 10.0, .v_max = v_max, .a_max = 10.0});
    auto ax2 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 20.0, .v_max = v_max, .a_max = 10.0});

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

    // The envelope is a hard limit, so the admissible slack is whatever rounding
    // the evaluation introduces -- and here that is none. The rebuilt profile has
    // a cruise plateau, the samples land on it, and on it the reported velocity
    // is the stored cruise velocity itself with only an exact unit sign factor
    // applied. So the sampled peak IS the profile's own reported peak, bitwise,
    // and that peak is at or below the limit the profile was rebuilt under. Two
    // exact comparisons replace a constant that was named for the machine epsilon
    // while being four hundred and fifty thousand times larger.
    double peak = 0.0;
    for (int i = 0; i <= 1000; ++i) {
        double const t = ax1.duration() * static_cast<double>(i) / 1000.0;
        auto const pt = ax1.evaluate(t);
        peak = std::max(peak, std::abs(pt.velocity[0]));
    }

    CAPTURE(peak, ax1.peak_velocity(), v_max);
    REQUIRE(peak == std::abs(ax1.peak_velocity()));
    REQUIRE(std::abs(ax1.peak_velocity()) <= v_max);
}

// --------------------------------------------------------------------------
// Test 6: Post-sync acceleration never exceeds a_max
// --------------------------------------------------------------------------
TEST_CASE("synchronize: post-sync acceleration within a_max", "[traj][sync]")
{
    double constexpr a_max = 10.0;
    auto ax1 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = a_max});
    auto ax2 = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = a_max});

    REQUIRE(ctrlpp::synchronize(ax1, ax2).has_value());

    // The acceleration needs no slack whatsoever, and the count that says so is
    // zero: the evaluation READS the stored acceleration magnitude, negates it on
    // the deceleration branch and applies a unit sign factor, so nothing along
    // that path can round. The retiming holds the magnitude fixed by
    // construction, and it is the commanded limit exactly for an axis commanded
    // from rest to rest. So the sampled peak is the limit itself, bitwise, and
    // never above it.
    double peak = 0.0;
    for (int i = 0; i <= 1000; ++i) {
        double const t = ax1.duration() * static_cast<double>(i) / 1000.0;
        auto const pt = ax1.evaluate(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= a_max);
        peak = std::max(peak, std::abs(pt.acceleration[0]));
    }

    CAPTURE(peak, a_max);
    REQUIRE(peak == a_max);
}

// --------------------------------------------------------------------------
// Test 7: Single axis is a no-op
// --------------------------------------------------------------------------
TEST_CASE("synchronize: single axis no-op", "[traj][sync]")
{
    auto ax = trapezoidal_axis(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    auto const before = snapshot(ax);

    REQUIRE(ctrlpp::synchronize(ax).has_value());

    // A single-axis synchronization is a no-op in the strongest sense: the target
    // is a bit-exact copy of this axis's own stored duration, so the retiming
    // compares equal and returns without touching anything. A tolerance there was
    // not merely unprincipled, it was the wrong kind of assertion -- it admitted a
    // rebuild that happened to land close. The file's own snapshot comparison
    // says what the case means, across the duration and nine interior samples of
    // position, velocity and acceleration, and the slowest-axis case at the end
    // already demonstrates that this exact form succeeds on this path.
    REQUIRE(identical(before, snapshot(ax)));
}

// --------------------------------------------------------------------------
// Test 8: Vector overload with trapezoidal axes
// --------------------------------------------------------------------------
TEST_CASE("synchronize: vector overload", "[traj][sync]")
{
    std::vector<ctrlpp::trapezoidal_trajectory<double>> axes;
    axes.push_back(trapezoidal_axis({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0}));
    axes.push_back(trapezoidal_axis({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0}));
    axes.push_back(trapezoidal_axis({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0}));

    double max_dur = 0.0;
    for (auto const& ax : axes) {
        max_dur = std::max(max_dur, ax.duration());
    }

    REQUIRE(ctrlpp::synchronize(std::span{axes}).has_value());

    // Same budget as the variadic overload: every axis here is trapezoidal, so
    // the trapezoidal rebuild count applies to each. The axis that already sets
    // the target lands on it bitwise, which the budget subsumes.
    double const budget = trapezoidal_duration_rounding_ops
                          * std::numeric_limits<double>::epsilon() * max_dur;
    for (auto const& ax : axes) {
        REQUIRE_THAT(ax.duration(), WithinAbs(max_dur, budget));
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
    auto faster = trapezoidal_axis(
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
