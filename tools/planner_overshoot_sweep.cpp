/// @file
/// @brief Probe: do the online planners' overshoot tests and the minimum
/// displacement a carry-velocity shape needs still disagree, and over what
/// population?
///
/// Two tests in the trajectory module answer overlapping questions about the
/// same command, with two different arithmetic chains behind them:
///
///   * the planner asks whether the distance it needs to stop exceeds the
///     distance that remains, and admits the command when it does not. The
///     comparison is relative between two lengths the planner has already
///     computed, with a slack that is a counted multiple of the scalar type's
///     epsilon;
///   * the shape the planner then asks for refuses a commanded displacement
///     below the distance the fastest admissible transition from the current
///     velocity to rest already sweeps. That refusal carries no slack at all.
///
/// Both quantities are the same length mathematically. They are not the same
/// number, because they are computed along different chains, and the planner's
/// side is additionally widened by its slack. Where the two boundaries differ
/// lies a band of commanded targets the planner admits and the shape cannot
/// serve. The planner reports every one of them, so the band is public behavior
/// rather than an internal detail.
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
/// **The two legs are not the same comparison, and the difference is a
/// finding.** The jerk-limited planner calls the double-S constructor directly,
/// so its leg compares the planner's verdict against that constructor's own
/// rejection on exactly the command the planner hands it. The trapezoidal
/// planner never calls that constructor -- it bounds no jerk and has no shape
/// with a domain to fall outside of -- so no such rejection exists to compare
/// against. Its leg compares the planner's overshoot verdict against the second
/// route the planner itself takes to the same length: the braking distance it
/// uses to place its own stopping point, formed as the mean speed times the
/// braking duration where the decision chain forms it as a squared speed over a
/// doubled acceleration limit. That is the closest analogue available, and it
/// is reported as an analogue rather than as the same measurement.
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
/// Building: standalone, no test framework and no build system. Eigen enters as
/// a system include, the way the library's own build treats it, so the only
/// diagnostics a build reports are this file's own.
///
///     g++ -std=c++20 -O2 -fno-exceptions -fno-rtti
///         -I lib/ctrlpp/include -isystem /usr/include/eigen3
///         tools/planner_overshoot_sweep.cpp -o /tmp/overshoot_sweep
///
/// Usage: no argument prints the full report. `--summary` prints the per
/// configuration measurements and the totals without the per-point lines.
///
/// **Where the boundaries this program locates were located.** A crossing that
/// moves with instruction selection is a crossing that has to name the
/// toolchain it was found on. This one does not move: the whole report is
/// BYTE-IDENTICAL across g++ 16.1.1, clang 22.1.8 and clang 18.1.8, at
/// optimization levels zero, two and three, and with floating-point contraction
/// off, at the compiler default and at `fast`, against Eigen 3.4.1. Eigen takes
/// no part in either verdict -- neither planner nor the double-S constructor
/// forms a matrix -- which is why one Eigen release is enough here where a
/// Riccati measurement would need two.

#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdint>
#include <cstring>

namespace {

using ctrlpp::double_s_trajectory;
using ctrlpp::online_planner_2nd;
using ctrlpp::online_planner_3rd;
using ctrlpp::online_planner_substitution_reason;
using ctrlpp::trajectory_error;

constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();

/// How far past each located boundary the exhaustive walk extends, in
/// representable targets. Wide enough that a boundary landing a few values away
/// from where the bisection put it is still inside the walk, and narrow enough
/// that the walk's own count stays readable.
constexpr std::int64_t walk_margin_ulps = 64;

/// The counted slacks the two planners carry, quoted so the measured relative
/// band width can be read against the width the counts predict. They enter no
/// comparison; nothing here is decided by them.
constexpr int third_order_slack_ops = 46;
constexpr int second_order_slack_ops = 4;

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
///
/// `slack_ops` is quoted beside the measured relative width so the two can be
/// read against each other. It enters no comparison.
template <typename PlannerVerdict, typename ShapeVerdict, typename Probe>
void measure_configuration(char const* configuration, double q0, double seed_distance,
                           double length_scale, int slack_ops, PlannerVerdict planner_verdict,
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
        std::printf("      band width relative   %.6e  against a counted slack of %.6e"
                    " (%d eps)\n",
                    (shape_side.target - planner_side.target) / remaining_at_planner,
                    static_cast<double>(slack_ops) * machine_epsilon, slack_ops);
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
                          third_order_slack_ops, planner_verdict, shape_verdict, probe,
                          print_points, leg);
}

// ---------------------------------------------------------------------------
// Second-order leg: the planner's overshoot test against the braking distance
// it computes for itself along a second route
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

    // The second route to the same length, and the planner's own: the distance
    // braking to rest sweeps, formed as the mean speed times the braking
    // duration, which is how the planner places its own stopping point. The
    // decision chain forms the same length as a squared speed over a doubled
    // acceleration limit. A commanded displacement below what braking already
    // sweeps has no carry-velocity trapezoid, for the same reason a commanded
    // displacement below the jerk-limited transition has no carry-velocity
    // double-S.
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
                          second_order_slack_ops, planner_verdict, shape_verdict, probe,
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
    {5.0, 10.0, 50.0},   // the set the coupled regression case commands
    {1.0, 5.0, 50.0},    // unit scale
    {0.25, 1.0, 4.0},    // slow axis
    {40.0, 20.0, 100.0}, // fast axis
    {1e-3, 1e-2, 1e-1},  // three decades down
    {1e3, 1e4, 1e5},     // three decades up
};

/// Axis positions at which the same limit set is probed. The band is a fixed
/// relative width, so its width in representable targets shrinks as the axis
/// travels away from the origin, and where it closes entirely a caller can no
/// longer land in it. That is measured rather than argued.
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

    std::printf("\n\nSecond-order planner: its overshoot test against the braking distance it\n");
    std::printf("computes for itself along a second route. The double-S constructor is not\n");
    std::printf("reachable from this planner, so no rejection of its own exists to compare\n");
    std::printf("against.\n\n");

    leg_counts second{};
    for (auto const& limits : second_order_limits) {
        for (auto const position : axis_positions) {
            for (auto const fraction : second_order_speed_fractions) {
                sweep_second_order_configuration(limits, position, fraction, print_points, second);
            }
        }
    }
    report_leg("Second-order", second);

    std::printf("\n");
    return 0;
}
