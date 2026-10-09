/// @file
/// @brief Probe: do the online planners' overshoot verdicts and the chain each
/// one reports against still disagree anywhere, and over what population? Does
/// every retargeted plan stay inside its limits, and does each shape
/// constructor honor its boundary-velocity precondition exactly?
///
/// The program runs two families of legs. The agreement legs, described first,
/// compare two booleans the library returns. The limit legs, described after
/// them, compare sampled values against the limits the caller stated.
///
/// Each planner decides whether a command overshoots, and then either carries
/// its current velocity into a profile or brakes to rest first. Each verdict is
/// meant to share one arithmetic chain with what it reports against:
///
///   * the jerk-limited planner's overshoot verdict calls the double-S shape's
///     own minimum-displacement predicate on the operands the shape's
///     constructor sees, so the planner and the constructor answer the same
///     question with the same arithmetic. A target the planner admits that the
///     shape still refuses is a disagreement, and the planner reports it as an
///     unavailable carry-velocity shape;
///   * the trapezoidal planner's overshoot verdict and the stopping point it
///     places are one braking distance, the mean speed times the braking
///     duration. This leg's reference is that documented formula, so it checks
///     the implementation against the formula rather than against a second
///     chain.
///
/// Where a verdict and its reference differ lies a band of commanded targets
/// the two classify differently. The planner reports every substitution, so a
/// band would be public behavior rather than an internal detail.
///
/// This program measures that band. It compares no quantity against a tolerance
/// of its own: both verdicts are booleans the library itself returns, and the
/// comparison between them is exact. A threshold inside the instrument would be
/// exactly the kind of underived constant this measurement exists to locate.
///
/// **What is measured, per configuration**
///
///   1. Each boundary is LOCATED rather than assumed, by bisecting the
///      commanded target in the space of representable doubles. It comes out as
///      the smallest representable target at which its verdict flips, so the
///      distance between the two is a count of representable targets -- units
///      in the last place of the commanded target -- and not a length that has
///      to be converted into one.
///   2. Every representable target from a fixed number of units in the last
///      place below the lower boundary to the same number above the upper one
///      is swept exhaustively. Nothing between the two is stepped over: the
///      band is walked one representable value at a time.
///   3. A geometric ladder of offsets, out to eleven decades of units in the
///      last place and eleven decades of the stopping distance in both
///      directions, establishes whether the disagreement is confined to the
///      band or is a property of the whole domain.
///
/// A verdict that does not flip anywhere in the widened bracket is REPORTED as
/// not flipping rather than dropped. A shape that refuses every commanded
/// displacement is a disagreement over the entire admissible domain, which is a
/// larger finding than a narrow band and would be invisible to an instrument
/// that required both boundaries to exist before it recorded anything.
///
/// The swept and disagreeing counts are exact numbers, and every disagreeing
/// point is printed with the configuration that produced it.
///
/// **The planner's verdict is read from the reason it reports, not from whether
/// it substituted.** A substitution reported as an unavailable carry-velocity
/// shape IS the shape's refusal showing through the planner, so reading the
/// disposition would compare the shape's verdict against itself and find
/// perfect agreement by construction.
///
/// **The two legs are not the same comparison.** The jerk-limited planner
/// calls the double-S constructor directly, so its leg compares the planner's
/// verdict against that constructor's own rejection on exactly the command the
/// planner hands it. The trapezoidal planner never calls that constructor -- it
/// bounds no jerk and has no shape with a domain to fall outside of -- so its
/// leg compares the planner's verdict against the braking distance formula the
/// planner documents for both its decision and its stopping point.
///
/// **Every state the sweep probes is one a caller can present.** The jerk
/// limited planner is driven to its own cruise phase and probed there, where
/// the sampled acceleration is exactly zero, so the acceleration-nulling phase
/// does not run and the command reaches the overshoot test carrying the
/// position and speed the planner itself reported. The shape constructor is
/// then handed those same two numbers. The trapezoidal planner carries no
/// acceleration state, so it is probed part way up its acceleration ramp and
/// its initial speed spans the reachable range. Neither leg reaches into a
/// private member or reconstructs one.
///
/// **The limit legs**
///
///   * Retarget envelope, in double and float. The jerk-limited planner is
///     driven from rest toward a far target and retargeted further along the
///     same direction from 40 instants in each of its cruise, acceleration-ramp
///     and acceleration-nulling windows, on eight limit sets. Every sample of
///     the new plan is checked against v_max and a_max: densely over its first
///     and last twenty seconds, where every phase but the cruise lives, and once
///     a second between. Every dense pair of samples is checked against j_max
///     times their spacing plus an allowance built from two counted rounding
///     chains, one on elapsed time and one on acceleration. That allowance is
///     printed beside j_max dt for each row, so the resolution of the jerk check
///     is visible: in float its time term dominates once elapsed time reaches
///     the order of a thousand seconds. A same-direction retarget that reports
///     any substitution is counted as a brake to rest that was not needed.
///   * Reversing, in double and float. The same planner is retargeted behind
///     its current position from the cruise and acceleration-ramp windows. It
///     must report a reversal and keep the same envelope while it brakes and
///     replans.
///   * Boundary-velocity precondition, for both the double-S and the
///     trapezoidal constructor in double and float. A boundary speed is
///     stepped from three representable values below v_max to six above it, as
///     the initial or the final velocity, in both directions of travel, under a
///     boundary-velocity rounding count of 0, 1, 2 and 4. A speed at or below
///     the limit must build. A speed above it must build only within the
///     admitted excess, and the clamped end must then report exactly +-v_max.
///     The realized peak, the double-S shape's reported peak or the trapezoidal
///     shape's peak sampled on a millisecond grid, must not exceed the limit.
///
/// Like the agreement legs, the limit legs reach no private member: they use
/// only create, update, sample, diagnostics, evaluate and peak_velocity.
///
/// Building: standalone, no test framework and no build system. Eigen enters as
/// a system include, the way the library's own build treats it, so the only
/// diagnostics a build reports are this file's own.
///
///     g++ -std=c++20 -O2 -Wall -Wextra -fno-exceptions -fno-rtti
///         -I lib/ctrlpp/include -isystem /usr/include/eigen3
///         tools/planner_overshoot_sweep.cpp -o build/tools/planner_overshoot_sweep
///
/// Usage: no argument prints the full report. `--summary` prints the per
/// configuration measurements and the totals without the per-point lines.
///
/// **Where the boundaries this program locates were located.** A crossing that
/// moves with instruction selection is a crossing that has to name the
/// toolchain it was found on. The archived report was produced by the build
/// command above with g++ 16.2.1 against Eigen 3.4.1, and by no other
/// toolchain. Eigen takes no part in either verdict -- neither planner nor the
/// double-S constructor forms a matrix.

#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <bit>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <type_traits>

namespace {

using ctrlpp::double_s_trajectory;
using ctrlpp::online_planner_2nd;
using ctrlpp::online_planner_3rd;
using ctrlpp::online_planner_substitution_reason;
using ctrlpp::trajectory_error;
using ctrlpp::trajectory_point;
using ctrlpp::trapezoidal_trajectory;

constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();

/// How far past each located boundary the exhaustive walk extends, in
/// representable targets. Wide enough that a boundary landing a few values away
/// from where the bisection put it is still inside the walk, and narrow enough
/// that the walk's own count stays readable.
constexpr std::int64_t walk_margin_ulps = 64;

// ---------------------------------------------------------------------------
// Arithmetic on representable targets
// ---------------------------------------------------------------------------

/// Position of a strictly positive double in the ordered sequence of
/// representable doubles. Adjacent representable values differ by one here, so
/// a difference of two such positions is a count of units in the last place
/// rather than a length that has to be divided by one.
auto target_index(double value) -> std::int64_t
{
    return static_cast<std::int64_t>(std::bit_cast<std::uint64_t>(value));
}

auto target_at(std::int64_t index) -> double
{
    return std::bit_cast<double>(static_cast<std::uint64_t>(index));
}

auto ulps_between(double low, double high) -> std::int64_t
{
    return target_index(high) - target_index(low);
}

// ---------------------------------------------------------------------------
// Names for what the library returned
// ---------------------------------------------------------------------------

auto error_name(trajectory_error error) -> char const*
{
    switch (error) {
    case trajectory_error::non_positive_velocity_limit:     return "non_positive_velocity_limit";
    case trajectory_error::non_positive_acceleration_limit: return "non_positive_acceleration_limit";
    case trajectory_error::non_positive_jerk_limit:         return "non_positive_jerk_limit";
    case trajectory_error::non_positive_duration:           return "non_positive_duration";
    case trajectory_error::non_finite_input:                return "non_finite_input";
    case trajectory_error::boundary_velocity_exceeds_limit: return "boundary_velocity_exceeds_limit";
    case trajectory_error::unrepresentable_duration:        return "unrepresentable_duration";
    case trajectory_error::unreachable_boundary_velocity:   return "unreachable_boundary_velocity";
    case trajectory_error::duration_shorter_than_current:   return "duration_shorter_than_current";
    case trajectory_error::unreachable_duration:            return "unreachable_duration";
    }
    return "unnamed";
}

auto reason_name(online_planner_substitution_reason reason) -> char const*
{
    switch (reason) {
    case online_planner_substitution_reason::none:                  return "none";
    case online_planner_substitution_reason::reversal_or_overshoot: return "reversal_or_overshoot";
    case online_planner_substitution_reason::carry_velocity_shape_unavailable:
        return "carry_velocity_shape_unavailable";
    }
    return "unnamed";
}

// ---------------------------------------------------------------------------
// What one probed command produced
// ---------------------------------------------------------------------------

/// The two verdicts at one commanded target, plus what the library reported
/// alongside them.
struct verdicts
{
    bool planner_overshoots{};  ///< the planner's own overshoot test selected the fallback
    bool shape_refuses{};       ///< the shape has no profile for this displacement
    online_planner_substitution_reason reason{online_planner_substitution_reason::none};
    trajectory_error shape_error{trajectory_error::non_finite_input};
    bool shape_error_valid{};
};

auto disagrees(verdicts const& result) -> bool
{
    return result.planner_overshoots != result.shape_refuses;
}

// ---------------------------------------------------------------------------
// Locating a boundary
// ---------------------------------------------------------------------------

/// Where a verdict turns over, or a statement that it does not.
///
/// The bracket is widened until the verdict holds at its low end and fails at
/// its high end, and both ends are then verified. A verdict that survives the
/// widening in either direction is recorded as such: that is a measurement, not
/// a failure of the search, and it turns out to be the shape of the largest
/// disagreement this sweep finds.
struct located_boundary
{
    double target{};
    bool located{};
    bool true_across_the_bracket{};
    bool false_across_the_bracket{};
    double bracket_low{};
    double bracket_high{};
};

template <typename Verdict>
auto locate_boundary(double q0, double seed_distance, Verdict verdict) -> located_boundary
{
    auto low_distance = seed_distance * 0.5;
    for (int attempt = 0; attempt < 60 && !verdict(q0 + low_distance); ++attempt) {
        low_distance *= 0.5;
    }
    auto high_distance = seed_distance * 4.0;
    for (int attempt = 0; attempt < 60 && verdict(q0 + high_distance); ++attempt) {
        high_distance *= 2.0;
    }

    located_boundary result{};
    result.bracket_low = q0 + low_distance;
    result.bracket_high = q0 + high_distance;

    if (!verdict(result.bracket_low)) {
        result.false_across_the_bracket = true;
        return result;
    }
    if (verdict(result.bracket_high)) {
        result.true_across_the_bracket = true;
        return result;
    }

    auto low = target_index(result.bracket_low);
    auto high = target_index(result.bracket_high);
    while (high - low > 1) {
        auto const middle = low + (high - low) / 2;
        if (verdict(target_at(middle))) {
            low = middle;
        } else {
            high = middle;
        }
    }
    result.target = target_at(high);
    result.located = true;
    return result;
}

// ---------------------------------------------------------------------------
// The exhaustive walk and the geometric ladder
// ---------------------------------------------------------------------------

/// What one configuration contributed.
struct configuration_counts
{
    long long walk_swept{};
    long long walk_disagreeing{};
    long long ladder_swept{};
    long long ladder_disagreeing{};
};

/// Powers-of-two offsets, in representable targets, at which the ladder probes
/// outside the walk window: from twice the walk margin out to eleven decades of
/// units in the last place away, in both directions.
constexpr int ladder_low_exponent = 7;
constexpr int ladder_high_exponent = 44;

/// Decades of the stopping distance the ladder additionally probes, so its
/// outer reaches are stated in the quantity the disagreement is about and not
/// in representable values alone.
constexpr int ladder_decades = 11;

template <typename Probe>
void walk_the_band(double lower, double upper, double reference, Probe probe,
                   configuration_counts& counts, bool print_points, char const* configuration)
{
    auto const first = target_index(lower) - walk_margin_ulps;
    auto const last = target_index(upper) + walk_margin_ulps;
    auto const origin = target_index(reference);

    for (auto index = first; index <= last; ++index) {
        auto const target = target_at(index);
        auto const result = probe(target);
        ++counts.walk_swept;
        if (!disagrees(result)) {
            continue;
        }
        ++counts.walk_disagreeing;
        if (!print_points) {
            continue;
        }
        std::printf("      DISAGREE  %s  walk %+7lld ulp  planner=%-10s shape=%-7s"
                    "  reported=%s  shape_error=%s\n",
                    configuration, static_cast<long long>(index - origin),
                    result.planner_overshoots ? "overshoots" : "admits",
                    result.shape_refuses ? "refuses" : "admits", reason_name(result.reason),
                    result.shape_error_valid ? error_name(result.shape_error) : "-");
    }
}

template <typename Probe>
void climb_the_ladder(double q0, double reference, double seed_distance, Probe probe,
                      configuration_counts& counts, bool print_points, char const* configuration)
{
    auto const origin = target_index(reference);

    // Only commands that still point at the target are recorded. A rung that
    // lands at or behind the start position asks a different question -- the
    // direction test's, not the overshoot test's -- and the two must not be
    // pooled into one population.
    auto const record = [&](double target, char const* rung) {
        if (!(target > q0)) {
            return;
        }
        auto const result = probe(target);
        ++counts.ladder_swept;
        if (!disagrees(result)) {
            return;
        }
        ++counts.ladder_disagreeing;
        if (!print_points) {
            return;
        }
        std::printf("      DISAGREE  %s  ladder %-6s %+lld ulp  planner=%-10s shape=%-7s"
                    "  reported=%s  shape_error=%s\n",
                    configuration, rung, static_cast<long long>(target_index(target) - origin),
                    result.planner_overshoots ? "overshoots" : "admits",
                    result.shape_refuses ? "refuses" : "admits", reason_name(result.reason),
                    result.shape_error_valid ? error_name(result.shape_error) : "-");
    };

    for (int exponent = ladder_low_exponent; exponent <= ladder_high_exponent; ++exponent) {
        auto const offset = std::int64_t{1} << exponent;
        record(target_at(origin + offset), "ulp");
        record(target_at(origin - offset), "ulp");
    }

    for (int decade = 0; decade <= ladder_decades; ++decade) {
        auto const distance = seed_distance * std::pow(10.0, -static_cast<double>(decade));
        record(reference + distance, "decade");
        record(reference - distance, "decade");
    }
}

// ---------------------------------------------------------------------------
// The measurement, shared by both legs
// ---------------------------------------------------------------------------

struct leg_counts
{
    long long walk_swept{};
    long long walk_disagreeing{};
    long long ladder_swept{};
    long long ladder_disagreeing{};
    int configurations{};
    int configurations_measured{};
    int configurations_with_a_band{};
    int configurations_refused_throughout{};
    std::int64_t widest_band_ulps{};
    std::int64_t narrowest_band_ulps{-1};
};

void accumulate(leg_counts& leg, configuration_counts const& counts)
{
    leg.walk_swept += counts.walk_swept;
    leg.walk_disagreeing += counts.walk_disagreeing;
    leg.ladder_swept += counts.ladder_swept;
    leg.ladder_disagreeing += counts.ladder_disagreeing;
}

/// Locate both boundaries, walk the region between them, climb the ladder
/// around them, and print what came out.
template <typename PlannerVerdict, typename ShapeVerdict, typename Probe>
void measure_configuration(char const* configuration, double q0, double seed_distance,
                           double length_scale, PlannerVerdict planner_verdict,
                           ShapeVerdict shape_verdict, Probe probe, bool print_points,
                           leg_counts& leg)
{
    auto const planner_side = locate_boundary(q0, seed_distance, planner_verdict);
    auto const shape_side = locate_boundary(q0, seed_distance, shape_verdict);

    if (!planner_side.located) {
        std::printf("      SKIPPED: the planner's overshoot verdict does not turn over across"
                    " remaining distances from %.6e to %.6e\n",
                    planner_side.bracket_low - q0, planner_side.bracket_high - q0);
        return;
    }

    auto const remaining_at_planner = planner_side.target - q0;
    std::printf("      planner boundary      %.17g  (%a)\n", planner_side.target,
                planner_side.target);
    std::printf("      remaining distance there   %.17g  = %.6e length scales\n",
                remaining_at_planner, remaining_at_planner / length_scale);

    configuration_counts counts{};
    ++leg.configurations_measured;

    if (shape_side.located) {
        auto const band_ulps = ulps_between(planner_side.target, shape_side.target);
        auto const magnitude = band_ulps < 0 ? -band_ulps : band_ulps;
        std::printf("      shape boundary        %.17g  (%a)\n", shape_side.target,
                    shape_side.target);
        std::printf("      band width            %lld ulp of the commanded target\n",
                    static_cast<long long>(band_ulps));
        std::printf("      band width relative   %.6e\n",
                    (shape_side.target - planner_side.target) / remaining_at_planner);
        if (band_ulps != 0) {
            ++leg.configurations_with_a_band;
        }
        if (magnitude > leg.widest_band_ulps) {
            leg.widest_band_ulps = magnitude;
        }
        if (leg.narrowest_band_ulps < 0 || magnitude < leg.narrowest_band_ulps) {
            leg.narrowest_band_ulps = magnitude;
        }
        if (band_ulps < 0) {
            std::printf("      NOTE: the shape's boundary sits BELOW the planner's, so the"
                        " disagreement runs the other way here\n");
            walk_the_band(shape_side.target, planner_side.target, planner_side.target, probe,
                          counts, print_points, configuration);
        } else {
            walk_the_band(planner_side.target, shape_side.target, planner_side.target, probe,
                          counts, print_points, configuration);
        }
    } else if (shape_side.true_across_the_bracket) {
        ++leg.configurations_refused_throughout;
        auto const at_high = probe(shape_side.bracket_high);
        std::printf("      shape boundary        NONE -- the shape refuses every commanded"
                    " displacement out to %.6e (%s)\n",
                    shape_side.bracket_high - q0,
                    at_high.shape_error_valid ? error_name(at_high.shape_error) : "-");
        std::printf("      band width            UNBOUNDED: every command the planner admits"
                    " is refused\n");
        walk_the_band(planner_side.target, planner_side.target, planner_side.target, probe, counts,
                      print_points, configuration);
    } else {
        std::printf("      shape boundary        NONE -- the shape admits every commanded"
                    " displacement down to %.6e\n",
                    shape_side.bracket_low - q0);
        walk_the_band(planner_side.target, planner_side.target, planner_side.target, probe, counts,
                      print_points, configuration);
    }

    climb_the_ladder(q0, planner_side.target, seed_distance, probe, counts, print_points,
                     configuration);

    std::printf("      swept %lld walk + %lld ladder = %lld;  disagreeing %lld walk + %lld"
                " ladder = %lld\n",
                counts.walk_swept, counts.ladder_swept, counts.walk_swept + counts.ladder_swept,
                counts.walk_disagreeing, counts.ladder_disagreeing,
                counts.walk_disagreeing + counts.ladder_disagreeing);
    accumulate(leg, counts);
}

// ---------------------------------------------------------------------------
// Third-order leg: the planner's overshoot test against the double-S
// constructor's own displacement rejection
// ---------------------------------------------------------------------------

struct limit_set
{
    double v_max{};
    double a_max{};
    double j_max{};
};

/// The state the planner reports part way through its own cruise phase.
///
/// Reported by the planner rather than assembled here. The acceleration is
/// required to be exactly zero, which is what keeps the acceleration-nulling
/// phase out of the command that follows: with it out, the command reaches the
/// overshoot test carrying the position and speed printed here, and the shape
/// constructor can be handed the same two numbers.
struct cruise_state
{
    double q0{};
    double v0{};
    double a0{};
    bool reached{};
};

auto drive_to_cruise(online_planner_3rd<double>& planner, limit_set const& limits,
                     double q_start) -> cruise_state
{
    planner.reset(q_start);

    auto const acceleration_stretch = limits.a_max / limits.j_max + limits.v_max / limits.a_max;
    auto const drive_time = 4.0 * acceleration_stretch;
    auto const far_target = q_start + limits.v_max * drive_time * 20.0;

    if (!planner.update(far_target).has_value()) {
        return {};
    }

    auto const point = planner.sample(drive_time);
    cruise_state state{};
    state.q0 = point.position[0];
    state.v0 = point.velocity[0];
    state.a0 = point.acceleration[0];
    state.reached = (state.a0 == 0.0) && (state.v0 > 0.0) && (state.q0 > 0.0);
    return state;
}

void sweep_third_order_configuration(limit_set const& limits, double q_start, bool print_points,
                                     leg_counts& leg)
{
    ++leg.configurations;

    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g j=%.4g q0~%.4g",
                  limits.v_max, limits.a_max, limits.j_max, q_start);

    auto made = online_planner_3rd<double>::create(
        {.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
    std::printf("  [%s]\n", configuration);
    if (!made.has_value()) {
        std::printf("      SKIPPED: the limits were rejected (%s)\n", error_name(made.error()));
        return;
    }
    auto& planner = *made;

    auto const cruise = drive_to_cruise(planner, limits, q_start);
    std::printf("      cruise state          q0=%.17g  v0=%.17g  a0=%.17g\n", cruise.q0, cruise.v0,
                cruise.a0);
    if (!cruise.reached) {
        std::printf("      SKIPPED: no cruise state with an exactly zero acceleration was"
                    " reached\n");
        return;
    }
    std::printf("      cruise speed against the velocity limit   %+lld ulp\n",
                static_cast<long long>(ulps_between(limits.v_max, cruise.v0)));
    std::printf("      deceleration reaches the acceleration limit   %s\n",
                (cruise.v0 > limits.a_max * limits.a_max / limits.j_max) ? "yes" : "no");

    // The per-point lines carry the speed the command was planned from as well
    // as the limits, so a disagreeing point is fully identified by its own line.
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g j=%.4g q0=%.17g v0=%.17g",
                  limits.v_max, limits.a_max, limits.j_max, cruise.q0, cruise.v0);

    // The planner's verdict is its own overshoot test, read from the reason it
    // reports. The shape's verdict is the constructor's, called on exactly the
    // configuration the planner hands it when it builds a carry-velocity
    // profile: the same start position, the same carried velocity, the same
    // limits and a terminal velocity of zero.
    auto const planner_probe = [&](double target) {
        verdicts result{};
        if (!planner.update(target).has_value()) {
            return result;
        }
        auto const& diagnostics = planner.diagnostics();
        result.reason = diagnostics.substitution_reason;
        result.planner_overshoots =
            (diagnostics.substitution_reason
             == online_planner_substitution_reason::reversal_or_overshoot);
        return result;
    };

    auto const shape_probe = [&](double target) {
        auto const built = double_s_trajectory<double>::create({.q0 = cruise.q0,
                                                               .q1 = target,
                                                               .v_max = limits.v_max,
                                                               .a_max = limits.a_max,
                                                               .j_max = limits.j_max,
                                                               .v0 = cruise.v0,
                                                               .v1 = 0.0});
        verdicts result{};
        if (built.has_value()) {
            return result;
        }
        result.shape_refuses = true;
        result.shape_error = built.error();
        result.shape_error_valid = true;
        return result;
    };

    auto const probe = [&](double target) {
        auto result = planner_probe(target);
        auto const shape = shape_probe(target);
        result.shape_refuses = shape.shape_refuses;
        result.shape_error = shape.shape_error;
        result.shape_error_valid = shape.shape_error_valid;
        return result;
    };

    auto const planner_verdict = [&](double target) {
        return planner_probe(target).planner_overshoots;
    };
    auto const shape_verdict = [&](double target) { return shape_probe(target).shape_refuses; };

    // A seed for the bracket only. It decides nothing: every bracket end is
    // widened until the verdict holds there and is then verified.
    auto const seed_distance = cruise.v0 * cruise.v0 / (2.0 * limits.a_max)
                               + limits.a_max * cruise.v0 / (2.0 * limits.j_max);

    // The direction test's length floor is a counted multiple of epsilon times
    // the stopping distance from full speed. Printing the remaining distance in
    // multiples of that scale shows the floor clause is not what selected the
    // substitution anywhere near the boundary.
    auto const length_scale = limits.v_max * limits.v_max / (2.0 * limits.a_max);

    measure_configuration(configuration, cruise.q0, seed_distance, length_scale,
                          planner_verdict, shape_verdict, probe,
                          print_points, leg);
}

// ---------------------------------------------------------------------------
// Second-order leg: the planner's overshoot test against the braking distance
// formula it documents for its decision and its stopping point
// ---------------------------------------------------------------------------

struct ramp_state
{
    double q0{};
    double v0{};
    bool reached{};
};

/// Sample the trapezoidal planner part way up its acceleration ramp.
///
/// This planner carries no acceleration state, so every speed on the ramp is a
/// state a caller can present to it, and the reachable speed range is swept
/// rather than probed at one point.
auto drive_to_ramp(online_planner_2nd<double>& planner, limit_set const& limits, double q_start,
                   double speed_fraction) -> ramp_state
{
    planner.reset(q_start);

    auto const ramp_time = limits.v_max / limits.a_max;
    auto const drive_time = speed_fraction * ramp_time;
    auto const far_target = q_start + limits.v_max * ramp_time * 200.0;

    if (!planner.update(far_target).has_value()) {
        return {};
    }

    auto const point = planner.sample(drive_time);
    ramp_state state{};
    state.q0 = point.position[0];
    state.v0 = point.velocity[0];
    state.reached = (state.v0 > 0.0) && (state.v0 <= limits.v_max) && (state.q0 > 0.0);
    return state;
}

void sweep_second_order_configuration(limit_set const& limits, double q_start,
                                      double speed_fraction, bool print_points, leg_counts& leg)
{
    ++leg.configurations;

    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g q0~%.4g v0/v=%.2f",
                  limits.v_max, limits.a_max, q_start, speed_fraction);

    auto made = online_planner_2nd<double>::create({.v_max = limits.v_max, .a_max = limits.a_max});
    std::printf("  [%s]\n", configuration);
    if (!made.has_value()) {
        std::printf("      SKIPPED: the limits were rejected (%s)\n", error_name(made.error()));
        return;
    }
    auto& planner = *made;

    auto const ramp = drive_to_ramp(planner, limits, q_start, speed_fraction);
    std::printf("      ramp state            q0=%.17g  v0=%.17g\n", ramp.q0, ramp.v0);
    if (!ramp.reached) {
        std::printf("      SKIPPED: no state on the acceleration ramp was reached\n");
        return;
    }

    // The per-point lines carry the speed the command was planned from as well
    // as the limits, so a disagreeing point is fully identified by its own line.
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g q0=%.17g v0=%.17g",
                  limits.v_max, limits.a_max, ramp.q0, ramp.v0);

    auto const planner_probe = [&](double target) {
        verdicts result{};
        if (!planner.update(target).has_value()) {
            return result;
        }
        auto const& diagnostics = planner.diagnostics();
        result.reason = diagnostics.substitution_reason;
        result.planner_overshoots =
            (diagnostics.substitution_reason
             == online_planner_substitution_reason::reversal_or_overshoot);
        return result;
    };

    // The distance braking to rest sweeps, formed as the mean speed times the
    // braking duration, which the planner documents as both its overshoot
    // verdict and its stopping point. A commanded displacement below it has no
    // carry-velocity trapezoid, for the same reason a commanded displacement
    // below the jerk-limited transition has no carry-velocity double-S.
    auto const shape_verdict = [&](double target) {
        auto const remaining = target - ramp.q0;
        auto const braking_duration = std::abs(ramp.v0) / limits.a_max;
        auto const braking_distance = std::abs(ramp.v0 * braking_duration / 2.0);
        return std::abs(remaining) < braking_distance;
    };

    auto const probe = [&](double target) {
        auto result = planner_probe(target);
        result.shape_refuses = shape_verdict(target);
        return result;
    };

    auto const planner_verdict = [&](double target) {
        return planner_probe(target).planner_overshoots;
    };

    auto const seed_distance = ramp.v0 * ramp.v0 / (2.0 * limits.a_max);
    auto const length_scale = limits.v_max * limits.v_max / (2.0 * limits.a_max);

    measure_configuration(configuration, ramp.q0, seed_distance, length_scale,
                          planner_verdict, shape_verdict, probe,
                          print_points, leg);
}

// ---------------------------------------------------------------------------
// The swept domain
// ---------------------------------------------------------------------------

/// Limit sets spanning several decades in every limit, chosen so that both
/// branches of the jerk-limited stopping-distance computation are covered: the
/// deceleration reaches the acceleration limit when the cruise speed exceeds
/// the squared acceleration limit over the jerk limit, and does not otherwise.
/// Which branch each set takes is printed with its result rather than asserted
/// here.
constexpr limit_set third_order_limits[] = {
    {5.0, 10.0, 50.0},   // the set the boundary-agreement regression case commands
    {1.0, 5.0, 50.0},    // unit scale
    {0.25, 1.0, 4.0},    // slow axis
    {40.0, 20.0, 100.0}, // fast axis
    {1e-3, 1e-2, 1e-1},  // three decades down
    {1e3, 1e4, 1e5},     // three decades up
};

/// Axis positions at which the same limit set is probed. A band of fixed
/// relative width shrinks in representable targets as the axis travels away
/// from the origin, so agreement is measured at several positions rather than
/// argued from one.
constexpr double axis_positions[] = {0.0, 1e2, 1e4};

constexpr limit_set second_order_limits[] = {
    {5.0, 10.0, 0.0}, {1.0, 5.0, 0.0}, {0.25, 1.0, 0.0}, {40.0, 20.0, 0.0}, {1e3, 1e4, 0.0},
};

constexpr double second_order_speed_fractions[] = {0.1, 0.35, 0.7, 0.95};

void report_leg(char const* name, leg_counts const& counts)
{
    std::printf("\n  %s totals\n", name);
    std::printf("      configurations                      %d\n", counts.configurations);
    std::printf("      configurations measured             %d\n", counts.configurations_measured);
    std::printf("      configurations with a bounded band  %d\n",
                counts.configurations_with_a_band);
    std::printf("      configurations refused throughout   %d\n",
                counts.configurations_refused_throughout);
    std::printf("      swept, exhaustive walk              %lld\n", counts.walk_swept);
    std::printf("      disagreeing, exhaustive walk        %lld\n", counts.walk_disagreeing);
    std::printf("      swept, geometric ladder             %lld\n", counts.ladder_swept);
    std::printf("      disagreeing, geometric ladder       %lld\n", counts.ladder_disagreeing);
    std::printf("      swept, total                        %lld\n",
                counts.walk_swept + counts.ladder_swept);
    std::printf("      disagreeing, total                  %lld\n",
                counts.walk_disagreeing + counts.ladder_disagreeing);
    std::printf("      widest bounded band                 %lld ulp\n",
                static_cast<long long>(counts.widest_band_ulps));
    std::printf("      narrowest bounded band              %lld ulp\n",
                static_cast<long long>(counts.narrowest_band_ulps));
}

// ---------------------------------------------------------------------------
// Retarget envelope and reversing legs: every sample of a retargeted plan
// against the limits the planner was given
// ---------------------------------------------------------------------------

// The last three are sets on which a cruise integrated forward through its ramp rounds above v_max.
constexpr limit_set retarget_limits[] = {
    {1.0, 1.0, 1.0},       {5.0, 10.0, 50.0}, {2.0, 3.0, 7.0},     {0.7, 1.3, 2.9},
    {10.0, 100.0, 1000.0}, {1.0, 5.0, 50.0},  {40.0, 20.0, 100.0}, {1e-3, 1e-2, 1e-1},
};

constexpr double grid_step = 1e-3;
constexpr std::int32_t instants_per_window = 40;
constexpr std::int32_t dense_seconds = 20;

// Elapsed time reaches a sample through at most 13 roundings, and the sampled acceleration carries at most 4
// since the last pinned phase end; a finite-difference jerk is admitted twice the excess those counts allow.
constexpr std::int32_t time_rounding_ops = 13;
constexpr std::int32_t acceleration_rounding_ops = 4;

template <typename Scalar>
struct typed_limits
{
    Scalar v_max;
    Scalar a_max;
    Scalar j_max;
};

template <typename Scalar>
auto typed(limit_set const& limits) -> typed_limits<Scalar>
{
    return {static_cast<Scalar>(limits.v_max), static_cast<Scalar>(limits.a_max),
            static_cast<Scalar>(limits.j_max)};
}

template <typename Scalar>
auto scalar_name() -> char const*
{
    return std::is_same_v<Scalar, float> ? "float" : "double";
}

enum class retarget_window : std::uint8_t
{
    cruise,
    acceleration_ramp,
    acceleration_nulling,
};

constexpr retarget_window retarget_windows[] = {retarget_window::cruise, retarget_window::acceleration_ramp,
                                                retarget_window::acceleration_nulling};

auto window_name(retarget_window window) -> char const*
{
    switch (window) {
    case retarget_window::cruise:               return "cruise";
    case retarget_window::acceleration_ramp:    return "acceleration ramp";
    case retarget_window::acceleration_nulling: return "acceleration nulling";
    }
    return "unnamed";
}

template <typename Scalar>
struct window_span
{
    Scalar begin;
    Scalar end;
};

// @cite biagiotti2009 -- Sec. 3.4.3, rest-to-rest ramp with and without a_max reached
template <typename Scalar>
auto window_span_of(retarget_window window, typed_limits<Scalar> const& limits, Scalar dt) -> window_span<Scalar>
{
    auto T_j = limits.a_max / limits.j_max;
    auto T_a = T_j + limits.v_max / limits.a_max;
    if (limits.v_max * limits.j_max < limits.a_max * limits.a_max) {
        T_j = std::sqrt(limits.v_max / limits.j_max);
        T_a = Scalar{2} * T_j;
    }
    switch (window) {
    case retarget_window::cruise:               return {Scalar{1.5} * T_a, Scalar{1.5} * T_a + Scalar{40} * dt};
    case retarget_window::acceleration_ramp:    return {dt, T_a - dt};
    case retarget_window::acceleration_nulling: break;
    }
    return {T_a - T_j, T_a - dt};
}

template <typename Scalar>
auto window_instant(window_span<Scalar> const& span, Scalar dt, std::int32_t i) -> Scalar
{
    auto const first = static_cast<std::int32_t>(std::ceil(span.begin / dt));
    auto const last = static_cast<std::int32_t>(std::floor(span.end / dt));
    auto const k = first + (last - first) * i / (instants_per_window - 1);
    return static_cast<Scalar>(k) * dt;
}

struct envelope_counts
{
    std::int64_t retargets;
    std::int64_t unexpected_reason;
    std::int64_t velocity_violations;
    std::int64_t acceleration_violations;
    std::int64_t jerk_violations;
};

void accumulate(envelope_counts& total, envelope_counts const& row)
{
    total.retargets += row.retargets;
    total.unexpected_reason += row.unexpected_reason;
    total.velocity_violations += row.velocity_violations;
    total.acceleration_violations += row.acceleration_violations;
    total.jerk_violations += row.jerk_violations;
}

template <typename Scalar>
struct envelope_scan
{
    online_planner_3rd<Scalar>& planner;
    typed_limits<Scalar> const& limits;
    envelope_counts& counts;
    char const* configuration;
    bool print_points;
    Scalar largest_allowance;
};

template <typename Scalar>
void report_violation(envelope_scan<Scalar> const& scan, char const* what, Scalar t, Scalar value)
{
    if (scan.print_points) {
        std::printf("      VIOLATION  [%s]  %s  t=%.17g  value=%.17g\n", scan.configuration, what,
                    static_cast<double>(t), static_cast<double>(value));
    }
}

template <typename Scalar>
void record_sample(envelope_scan<Scalar>& scan, Scalar t, trajectory_point<Scalar, 1> const& point)
{
    auto const v = std::abs(point.velocity[0]);
    auto const a = std::abs(point.acceleration[0]);
    if (v > scan.limits.v_max) {
        ++scan.counts.velocity_violations;
        report_violation(scan, "|v| above v_max", t, v);
    }
    if (a > scan.limits.a_max) {
        ++scan.counts.acceleration_violations;
        report_violation(scan, "|a| above a_max", t, a);
    }
}

template <typename Scalar>
void record_jerk(envelope_scan<Scalar>& scan, Scalar t_prev, Scalar t, Scalar a_prev, Scalar a)
{
    auto const eps = std::numeric_limits<Scalar>::epsilon();
    auto const& limits = scan.limits;
    auto const time_scale = std::max(std::abs(t), scan.planner.diagnostics().planned_duration);
    auto const allowance = Scalar{2}
                           * (static_cast<Scalar>(time_rounding_ops) * eps * limits.j_max * time_scale
                              + static_cast<Scalar>(acceleration_rounding_ops) * eps * limits.a_max);
    scan.largest_allowance = std::max(scan.largest_allowance, allowance);
    auto const bound = limits.j_max * (t - t_prev) + allowance;
    if (std::abs(a - a_prev) > bound) {
        ++scan.counts.jerk_violations;
        report_violation(scan, "|a_k - a_(k-1)| above j_max dt plus the allowance", t, std::abs(a - a_prev) - bound);
    }
}

template <typename Scalar>
void scan_dense(envelope_scan<Scalar>& scan, Scalar t_begin)
{
    auto const dt = static_cast<Scalar>(grid_step);
    auto const steps = static_cast<std::int32_t>(std::lround(dense_seconds / grid_step));
    Scalar t_prev{};
    Scalar a_prev{};
    for (std::int32_t k = 0; k <= steps; ++k) {
        auto const t = t_begin + static_cast<Scalar>(k) * dt;
        auto const point = scan.planner.sample(t);
        record_sample(scan, t, point);
        if (k > 0) {
            record_jerk(scan, t_prev, t, a_prev, point.acceleration[0]);
        }
        t_prev = t;
        a_prev = point.acceleration[0];
    }
}

// Dense over the first and last stretch of the plan, where every phase but the cruise lives, and every second between.
template <typename Scalar>
void scan_plan(envelope_scan<Scalar>& scan, Scalar t0)
{
    auto const duration = scan.planner.diagnostics().planned_duration;
    auto const dense = static_cast<Scalar>(dense_seconds);
    scan_dense(scan, t0);
    for (std::int32_t k = dense_seconds; static_cast<Scalar>(k) + dense < duration; ++k) {
        auto const t = t0 + static_cast<Scalar>(k);
        record_sample(scan, t, scan.planner.sample(t));
    }
    scan_dense(scan, t0 + std::max(duration - dense, Scalar{0}));
}

template <typename Scalar>
void record_reason(envelope_scan<Scalar>& scan, Scalar t, bool reversing)
{
    auto const expected = reversing ? online_planner_substitution_reason::reversal_or_overshoot
                                    : online_planner_substitution_reason::none;
    auto const reason = scan.planner.diagnostics().substitution_reason;
    if (reason == expected) {
        return;
    }
    ++scan.counts.unexpected_reason;
    if (scan.print_points) {
        std::printf("      REASON  [%s]  t=%.17g  reported=%s\n", scan.configuration, static_cast<double>(t),
                    reason_name(reason));
    }
}

template <typename Scalar>
void retarget_once(envelope_scan<Scalar>& scan, Scalar t, bool reversing)
{
    auto& planner = scan.planner;
    auto const v_max = scan.limits.v_max;
    planner.reset(Scalar{0});
    auto const driven = planner.update(Scalar{1000} * v_max).has_value();
    auto const q = planner.sample(t).position[0];
    auto const target = reversing ? q - Scalar{1000} * v_max : Scalar{2000} * v_max;
    if (!driven || !planner.update(target).has_value()) {
        std::printf("      REFUSED  [%s]  t=%.17g\n", scan.configuration, static_cast<double>(t));
        return;
    }
    ++scan.counts.retargets;
    record_reason(scan, t, reversing);
    scan_plan(scan, t);
}

template <typename Scalar>
void measure_window(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits,
                    retarget_window window, bool reversing, bool print_points, envelope_counts& total)
{
    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g j=%.4g %s", static_cast<double>(limits.v_max),
                  static_cast<double>(limits.a_max), static_cast<double>(limits.j_max), window_name(window));
    envelope_counts row{};
    envelope_scan<Scalar> scan{planner, limits, row, configuration, print_points, Scalar{0}};
    auto const dt = static_cast<Scalar>(grid_step);
    auto const span = window_span_of(window, limits, dt);
    for (std::int32_t i = 0; i < instants_per_window; ++i) {
        retarget_once(scan, window_instant(span, dt, i), reversing);
    }
    std::printf("  [%s]  retargets %lld  reason %lld  |v|>v_max %lld  |a|>a_max %lld  jerk %lld"
                "  allowance %.3e  j_max dt %.3e\n",
                configuration, static_cast<long long>(row.retargets), static_cast<long long>(row.unexpected_reason),
                static_cast<long long>(row.velocity_violations), static_cast<long long>(row.acceleration_violations),
                static_cast<long long>(row.jerk_violations), static_cast<double>(scan.largest_allowance),
                static_cast<double>(limits.j_max * dt));
    accumulate(total, row);
}

void report_count(char const* label, std::int64_t count)
{
    std::printf("      %-44s %lld\n", label, static_cast<long long>(count));
}

void report_envelope(char const* leg, char const* type, char const* reason_label, envelope_counts const& counts)
{
    std::printf("\n  %s, %s totals\n", leg, type);
    report_count("retargets", counts.retargets);
    report_count(reason_label, counts.unexpected_reason);
    report_count("velocity samples above the limit", counts.velocity_violations);
    report_count("acceleration samples above the limit", counts.acceleration_violations);
    report_count("jerk pairs above the allowance", counts.jerk_violations);
}

template <typename Scalar>
void run_retarget_leg(char const* leg, bool reversing, bool print_points)
{
    std::printf("\n%s, %s\n\n", leg, scalar_name<Scalar>());
    envelope_counts total{};
    auto const windows = reversing ? 2 : 3;
    for (auto const& raw : retarget_limits) {
        auto const limits = typed<Scalar>(raw);
        auto made = online_planner_3rd<Scalar>::create(
            {.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
        if (!made.has_value()) {
            std::printf("  [v=%.4g]  SKIPPED: %s\n", raw.v_max, error_name(made.error()));
            continue;
        }
        for (auto w = 0; w < windows; ++w) {
            measure_window(*made, limits, retarget_windows[w], reversing, print_points, total);
        }
    }
    report_envelope(leg, scalar_name<Scalar>(), reversing ? "reason other than reversal_or_overshoot" : "substituted",
                    total);
}

// ---------------------------------------------------------------------------
// Boundary-velocity precondition of both shape constructors
// ---------------------------------------------------------------------------

constexpr std::int32_t lowest_speed_step = -3;
constexpr std::int32_t highest_speed_step = 6;
constexpr std::uint32_t knob_values[] = {0, 1, 2, 4};
constexpr auto error_count = static_cast<std::size_t>(trajectory_error::unreachable_duration) + 1;

struct precondition_counts
{
    std::int64_t commands;
    std::int64_t builds;
    std::array<std::int64_t, error_count> refusals;
    std::int64_t peak_above_limit;
    std::int64_t refused_at_or_below;
    std::int64_t refused_within_allowance;
    std::int64_t admitted_beyond_allowance;
    std::int64_t clamped_exact;
    std::int64_t clamped_off;
    double largest_peak_ratio;
};

void accumulate(precondition_counts& total, precondition_counts const& row)
{
    total.commands += row.commands;
    total.builds += row.builds;
    for (std::size_t e = 0; e < error_count; ++e) {
        total.refusals[e] += row.refusals[e];
    }
    total.peak_above_limit += row.peak_above_limit;
    total.refused_at_or_below += row.refused_at_or_below;
    total.refused_within_allowance += row.refused_within_allowance;
    total.admitted_beyond_allowance += row.admitted_beyond_allowance;
    total.clamped_exact += row.clamped_exact;
    total.clamped_off += row.clamped_off;
    total.largest_peak_ratio = std::max(total.largest_peak_ratio, row.largest_peak_ratio);
}

template <typename Scalar>
struct boundary_command
{
    std::int32_t step;
    Scalar speed;
    Scalar sigma;
    bool final_role;
    std::uint32_t knob;
};

template <typename Scalar>
auto stepped_speed(Scalar v_max, std::int32_t step) -> Scalar
{
    auto const toward = step < 0 ? Scalar{0} : std::numeric_limits<Scalar>::infinity();
    auto speed = v_max;
    for (std::int32_t k = 0; k < std::abs(step); ++k) {
        speed = std::nextafter(speed, toward);
    }
    return speed;
}

template <typename Scalar>
auto allowance_of(boundary_command<Scalar> const& command, Scalar v_max) -> Scalar
{
    return static_cast<Scalar>(command.knob) * std::numeric_limits<Scalar>::epsilon() * v_max;
}

template <typename Shape, typename Scalar>
auto build_shape(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command)
    -> ctrlpp::expected<Shape, trajectory_error>
{
    auto const v = command.sigma * command.speed;
    auto const q1 = command.sigma * Scalar{10} * limits.v_max;
    auto const v0 = command.final_role ? Scalar{0} : v;
    auto const v1 = command.final_role ? v : Scalar{0};
    if constexpr (std::is_same_v<Shape, double_s_trajectory<Scalar>>) {
        return Shape::create({.q0 = Scalar{0}, .q1 = q1, .v_max = limits.v_max, .a_max = limits.a_max,
                              .j_max = limits.j_max, .v0 = v0, .v1 = v1,
                              .boundary_velocity_rounding_ops = command.knob});
    } else {
        return Shape::create({.q0 = Scalar{0}, .q1 = q1, .v_max = limits.v_max, .a_max = limits.a_max,
                              .v0 = v0, .v1 = v1, .boundary_velocity_rounding_ops = command.knob});
    }
}

template <typename Scalar>
auto realized_peak(double_s_trajectory<Scalar> const& shape) -> Scalar
{
    return std::abs(shape.peak_velocity());
}

template <typename Scalar>
auto realized_peak(trapezoidal_trajectory<Scalar> const& shape) -> Scalar
{
    auto const dt = static_cast<Scalar>(grid_step);
    auto peak = std::abs(shape.evaluate(shape.duration()).velocity[0]);
    for (std::int32_t k = 0; static_cast<Scalar>(k) * dt < shape.duration(); ++k) {
        peak = std::max(peak, std::abs(shape.evaluate(static_cast<Scalar>(k) * dt).velocity[0]));
    }
    return peak;
}

template <typename Scalar>
void print_command(char const* configuration, boundary_command<Scalar> const& command, char const* outcome)
{
    std::printf("      [%s  speed %+d ulp  sigma %+d  %s  knob %u]  %s\n", configuration, command.step,
                static_cast<std::int32_t>(command.sigma), command.final_role ? "v1" : "v0", command.knob, outcome);
}

template <typename Scalar>
void record_refusal(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command,
                    trajectory_error error, precondition_counts& counts)
{
    ++counts.refusals[static_cast<std::size_t>(error)];
    if (command.speed <= limits.v_max) {
        ++counts.refused_at_or_below;
    } else if (command.speed - limits.v_max <= allowance_of(command, limits.v_max)) {
        ++counts.refused_within_allowance;
    }
}

template <typename Shape, typename Scalar>
auto record_build(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command, Shape const& shape,
                  precondition_counts& counts) -> Scalar
{
    ++counts.builds;
    auto const peak = realized_peak(shape);
    counts.peak_above_limit += (peak > limits.v_max) ? 1 : 0;
    counts.largest_peak_ratio = std::max(counts.largest_peak_ratio, static_cast<double>(peak / limits.v_max));
    if (command.speed > limits.v_max) {
        counts.admitted_beyond_allowance += (command.speed - limits.v_max > allowance_of(command, limits.v_max)) ? 1 : 0;
        auto const end = shape.evaluate(command.final_role ? shape.duration() : Scalar{0}).velocity[0];
        ++(end == command.sigma * limits.v_max ? counts.clamped_exact : counts.clamped_off);
    }
    return peak;
}

template <typename Shape, typename Scalar>
void measure_command(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command,
                     char const* configuration, bool print_points, precondition_counts& counts)
{
    ++counts.commands;
    char outcome[96]{};
    auto const built = build_shape<Shape>(limits, command);
    if (built.has_value()) {
        auto const peak = record_build(limits, command, *built, counts);
        std::snprintf(outcome, sizeof(outcome), "built  peak %.17g", static_cast<double>(peak));
    } else {
        record_refusal(limits, command, built.error(), counts);
        std::snprintf(outcome, sizeof(outcome), "refused  %s", error_name(built.error()));
    }
    if (print_points) {
        print_command(configuration, command, outcome);
    }
}

void print_precondition_row(char const* configuration, precondition_counts const& row)
{
    std::printf("  [%s]  commands %lld  builds %lld", configuration, static_cast<long long>(row.commands),
                static_cast<long long>(row.builds));
    for (std::size_t e = 0; e < error_count; ++e) {
        if (row.refusals[e] != 0) {
            std::printf("  %s %lld", error_name(static_cast<trajectory_error>(e)),
                        static_cast<long long>(row.refusals[e]));
        }
    }
    std::printf("  largest peak / v_max - 1 %+.3e\n", row.largest_peak_ratio - 1.0);
}

template <typename Shape, typename Scalar>
void measure_limit_set(typed_limits<Scalar> const& limits, char const* shape_name, bool print_points,
                       precondition_counts& total)
{
    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "%s %s v=%.4g a=%.4g j=%.4g", shape_name,
                  scalar_name<Scalar>(), static_cast<double>(limits.v_max), static_cast<double>(limits.a_max),
                  static_cast<double>(limits.j_max));
    precondition_counts row{};
    for (auto step = lowest_speed_step; step <= highest_speed_step; ++step) {
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            for (auto const final_role : {false, true}) {
                for (auto const knob : knob_values) {
                    boundary_command<Scalar> const command{step, stepped_speed(limits.v_max, step), sigma,
                                                           final_role, knob};
                    measure_command<Shape>(limits, command, configuration, print_points, row);
                }
            }
        }
    }
    print_precondition_row(configuration, row);
    accumulate(total, row);
}

void report_precondition(char const* shape_name, char const* type, precondition_counts const& counts)
{
    std::printf("\n  Boundary-velocity precondition, %s, %s totals\n", shape_name, type);
    report_count("commands", counts.commands);
    report_count("built", counts.builds);
    for (std::size_t e = 0; e < error_count; ++e) {
        if (counts.refusals[e] != 0) {
            std::printf("      refused, %-35s %lld\n", error_name(static_cast<trajectory_error>(e)),
                        static_cast<long long>(counts.refusals[e]));
        }
    }
    report_count("builds with a peak above the limit", counts.peak_above_limit);
    report_count("refusals at or below the limit", counts.refused_at_or_below);
    report_count("refused within the allowance above the limit", counts.refused_within_allowance);
    report_count("admissions beyond the allowance", counts.admitted_beyond_allowance);
    report_count("clamped ends reported exactly at the limit", counts.clamped_exact);
    report_count("clamped ends reported off the limit", counts.clamped_off);
}

template <typename Shape>
void run_precondition_leg(char const* shape_name, bool print_points)
{
    using Scalar = typename Shape::scalar_type;
    std::printf("\nBoundary-velocity precondition, %s, %s\n\n", shape_name, scalar_name<Scalar>());
    precondition_counts total{};
    for (auto const& raw : retarget_limits) {
        measure_limit_set<Shape>(typed<Scalar>(raw), shape_name, print_points, total);
    }
    report_precondition(shape_name, scalar_name<Scalar>(), total);
}

void run_limit_legs(bool print_points)
{
    std::printf("\n\nThird-order planner retargeted from its cruise, acceleration-ramp and acceleration-nulling\n");
    std::printf("windows: every sample of the new plan against v_max and a_max, and every dense pair of\n");
    std::printf("samples against j_max plus the counted rounding allowance.\n");
    run_retarget_leg<double>("Retarget envelope", false, print_points);
    run_retarget_leg<float>("Retarget envelope", false, print_points);
    run_retarget_leg<double>("Reversing", true, print_points);
    run_retarget_leg<float>("Reversing", true, print_points);

    std::printf("\n\nShape constructors handed a boundary speed stepped across v_max one representable value\n");
    std::printf("at a time, as v0 and as v1, in both frames, under each boundary-velocity rounding count.\n");
    run_precondition_leg<double_s_trajectory<double>>("double_s_trajectory", print_points);
    run_precondition_leg<double_s_trajectory<float>>("double_s_trajectory", print_points);
    run_precondition_leg<trapezoidal_trajectory<double>>("trapezoidal_trajectory", print_points);
    run_precondition_leg<trapezoidal_trajectory<float>>("trapezoidal_trajectory", print_points);
}

}

auto main(int argc, char** argv) -> int
{
    bool print_points = true;
    for (int index = 1; index < argc; ++index) {
        if (std::strcmp(argv[index], "--summary") == 0) {
            print_points = false;
        }
    }

    std::printf("Planner overshoot verdict against the minimum displacement a carry-velocity"
                " shape needs\n");
    std::printf("double, epsilon = %.6e\n\n", machine_epsilon);

    std::printf("Third-order planner: its overshoot test against the double-S constructor's own\n");
    std::printf("displacement rejection, on exactly the command the planner hands it.\n\n");

    leg_counts third{};
    for (auto const& limits : third_order_limits) {
        for (auto const position : axis_positions) {
            sweep_third_order_configuration(limits, position, print_points, third);
        }
    }
    report_leg("Third-order", third);

    std::printf("\n\nSecond-order planner: its overshoot test against the braking distance\n");
    std::printf("formula it documents for its decision and its stopping point. The double-S\n");
    std::printf("constructor is not reachable from this planner, so no rejection of its own\n");
    std::printf("exists to compare against.\n\n");

    leg_counts second{};
    for (auto const& limits : second_order_limits) {
        for (auto const position : axis_positions) {
            for (auto const fraction : second_order_speed_fractions) {
                sweep_second_order_configuration(limits, position, fraction, print_points, second);
            }
        }
    }
    report_leg("Second-order", second);

    run_limit_legs(print_points);

    std::printf("\n");
    return 0;
}
