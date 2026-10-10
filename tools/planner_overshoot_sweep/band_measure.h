#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_BAND_MEASURE_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_BAND_MEASURE_H

#include "common.h"
#include "boundary.h"
#include "band_walk.h"

#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cstdio>
#include <cstdint>
#include <algorithm>

namespace planner_overshoot_sweep
{

struct leg_counts
{
    long long walk_swept{};
    long long walk_disagreeing{};
    long long ladder_swept{};
    long long ladder_disagreeing{};
    long long walk_behind_start{};
    long long behind_start_not_reversal{};
    int configurations{};
    int configurations_measured{};
    int configurations_with_a_band{};
    int configurations_refused_throughout{};
    std::int64_t widest_band_ulps{};
    std::int64_t narrowest_band_ulps{-1};
};

inline void accumulate(leg_counts& leg, configuration_counts const& counts)
{
    leg.walk_swept += counts.walk_swept;
    leg.walk_disagreeing += counts.walk_disagreeing;
    leg.ladder_swept += counts.ladder_swept;
    leg.ladder_disagreeing += counts.ladder_disagreeing;
    leg.walk_behind_start += counts.walk_behind_start;
    leg.behind_start_not_reversal += counts.behind_start_not_reversal;
}

inline void record_band(leg_counts& leg, std::int64_t band_ulps)
{
    auto const magnitude = band_ulps < 0 ? -band_ulps : band_ulps;
    if (band_ulps != 0) {
        ++leg.configurations_with_a_band;
    }
    leg.widest_band_ulps = std::max(leg.widest_band_ulps, magnitude);
    if (leg.narrowest_band_ulps < 0 || magnitude < leg.narrowest_band_ulps) {
        leg.narrowest_band_ulps = magnitude;
    }
}

template <typename Scalar>
void print_planner_boundary(configuration_frame<Scalar> const& frame, Scalar target)
{
    auto const remaining = target - frame.q0;
    std::printf("      planner boundary      %.17g  (%a)\n", static_cast<double>(target), static_cast<double>(target));
    std::printf("      remaining distance there   %.17g  = %.6e length scales\n", static_cast<double>(remaining),
                static_cast<double>(remaining / frame.length_scale));
}

template <typename Scalar, typename Probe>
void walk_bounded_band(configuration_frame<Scalar> const& frame, Scalar planner_target, Scalar shape_target,
                       Probe probe, configuration_counts& counts, leg_counts& leg)
{
    auto const band_ulps = static_cast<std::int64_t>(ulps_between(planner_target, shape_target));
    std::printf("      shape boundary        %.17g  (%a)\n", static_cast<double>(shape_target),
                static_cast<double>(shape_target));
    std::printf("      band width            %lld ulp of the commanded target\n", static_cast<long long>(band_ulps));
    std::printf("      band width relative   %.6e\n",
                static_cast<double>((shape_target - planner_target) / (planner_target - frame.q0)));
    record_band(leg, band_ulps);
    if (band_ulps < 0) {
        std::printf("      NOTE: the shape's boundary sits BELOW the planner's, so the"
                    " disagreement runs the other way here\n");
        walk_the_band(frame, shape_target, planner_target, planner_target, probe, counts);
    } else {
        walk_the_band(frame, planner_target, shape_target, planner_target, probe, counts);
    }
}

template <typename Scalar, typename Probe>
void walk_unbounded_band(configuration_frame<Scalar> const& frame, located_boundary<Scalar> const& shape_side,
                         Scalar planner_target, Probe probe, configuration_counts& counts, leg_counts& leg)
{
    if (shape_side.true_across_the_bracket) {
        ++leg.configurations_refused_throughout;
        auto const at_high = probe(shape_side.bracket_high);
        std::printf("      shape boundary        NONE -- the shape refuses every commanded"
                    " displacement out to %.6e (%s)\n",
                    static_cast<double>(shape_side.bracket_high - frame.q0),
                    at_high.shape_error_valid ? error_name(at_high.shape_error) : "-");
        std::printf("      band width            UNBOUNDED: every command the planner admits"
                    " is refused\n");
    } else {
        std::printf("      shape boundary        NONE -- the shape admits every commanded"
                    " displacement down to %.6e\n",
                    static_cast<double>(shape_side.bracket_low - frame.q0));
    }
    walk_the_band(frame, planner_target, planner_target, planner_target, probe, counts);
}

inline void print_configuration_counts(configuration_counts const& counts)
{
    std::printf("      swept %lld walk + %lld ladder = %lld;  disagreeing %lld walk + %lld"
                " ladder = %lld\n",
                counts.walk_swept, counts.ladder_swept, counts.walk_swept + counts.ladder_swept,
                counts.walk_disagreeing, counts.ladder_disagreeing,
                counts.walk_disagreeing + counts.ladder_disagreeing);
    std::printf("      behind the start      %lld walk targets skipped;  not reported as reversal_or_overshoot %lld\n",
                counts.walk_behind_start, counts.behind_start_not_reversal);
}

template <typename Scalar>
void report_unturned(configuration_frame<Scalar> const& frame, located_boundary<Scalar> const& planner_side)
{
    std::printf("      SKIPPED: the planner's overshoot verdict does not turn over across"
                " remaining distances from %.6e to %.6e\n",
                static_cast<double>(planner_side.bracket_low - frame.q0),
                static_cast<double>(planner_side.bracket_high - frame.q0));
}

/// Locate both boundaries, walk the region between them, climb the ladder
/// around them, and print what came out.
template <typename Scalar, typename PlannerVerdict, typename ShapeVerdict, typename Probe>
void measure_configuration(configuration_frame<Scalar> const& frame, PlannerVerdict planner_verdict,
                           ShapeVerdict shape_verdict, Probe probe, leg_counts& leg)
{
    auto const planner_side = locate_boundary(frame.q0, frame.seed_distance, planner_verdict);
    auto const shape_side = locate_boundary(frame.q0, frame.seed_distance, shape_verdict);
    if (!planner_side.located) {
        report_unturned(frame, planner_side);
        return;
    }
    print_planner_boundary(frame, planner_side.target);
    configuration_counts counts{};
    ++leg.configurations_measured;
    if (shape_side.located) {
        walk_bounded_band(frame, planner_side.target, shape_side.target, probe, counts, leg);
    } else {
        walk_unbounded_band(frame, shape_side, planner_side.target, probe, counts, leg);
    }
    climb_the_ladder(frame.q0, planner_side.target, frame.seed_distance, probe, counts, frame.print_points,
                     frame.configuration);
    print_configuration_counts(counts);
    accumulate(leg, counts);
}

/// The planner's verdict is its own overshoot test, read from the reason it
/// reports.
template <typename Planner, typename Scalar>
auto planner_verdicts(Planner& planner, Scalar target) -> verdicts
{
    verdicts result{};
    if (!planner.update(target).has_value()) {
        return result;
    }
    auto const& diagnostics = planner.diagnostics();
    result.reason = diagnostics.substitution_reason;
    result.planner_overshoots =
        (diagnostics.substitution_reason == online_planner_substitution_reason::reversal_or_overshoot);
    return result;
}

}

#endif
