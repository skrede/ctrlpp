#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_BAND_WALK_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_BAND_WALK_H

#include "common.h"
#include "boundary.h"

#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdint>

namespace planner_overshoot_sweep
{

/// What one configuration contributed.
struct configuration_counts
{
    long long walk_swept{};
    long long walk_disagreeing{};
    long long ladder_swept{};
    long long ladder_disagreeing{};
    long long walk_behind_start{};
    long long behind_start_not_reversal{};
};

/// The command a configuration was planned from, and how its report is printed.
template <typename Scalar>
struct configuration_frame
{
    char const* configuration;
    Scalar q0;
    Scalar seed_distance;
    Scalar length_scale;
    bool print_points;
};

/// Powers-of-two offsets, in representable targets, at which the ladder probes
/// outside the walk window, in both directions, from twice the walk margin out.
inline constexpr std::int32_t ladder_low_exponent = 7;

// 2^(digits - 9) units in the last place move a target by 2^-9 to 2^-8 of itself, inside the significand of either type.
template <typename Scalar>
constexpr std::int32_t ladder_high_exponent = std::numeric_limits<Scalar>::digits - 9;

/// Decades of the stopping distance the ladder additionally probes, so its
/// outer reaches are stated in the quantity the disagreement is about and not
/// in representable values alone.
inline constexpr std::int32_t ladder_decades = 11;

inline void record_walk(verdicts const& result, char const* where, configuration_counts& counts, bool print_points,
                        char const* configuration)
{
    ++counts.walk_swept;
    if (!disagrees(result)) {
        return;
    }
    ++counts.walk_disagreeing;
    if (print_points) {
        print_disagreement(configuration, where, result);
    }
}

// A command at or behind the start position is the direction test's question, so it leaves the overshoot pool and
// is held to the direction test's answer instead.
inline void record_behind_start(verdicts const& result, char const* where, configuration_counts& counts, bool print_points,
                                char const* configuration)
{
    ++counts.walk_behind_start;
    if (result.reason == online_planner_substitution_reason::reversal_or_overshoot) {
        return;
    }
    ++counts.behind_start_not_reversal;
    if (print_points) {
        std::printf("      BEHIND-START  %s  %s  reported=%s\n", configuration, where, reason_name(result.reason));
    }
}

template <typename Scalar, typename Probe>
void walk_the_band(configuration_frame<Scalar> const& frame, Scalar lower, Scalar upper, Scalar reference,
                   Probe probe, configuration_counts& counts)
{
    auto const first = target_index(lower) - walk_margin_ulps;
    auto const last = target_index(upper) + walk_margin_ulps;
    auto const origin = target_index(reference);
    for (auto index = first; index <= last; ++index) {
        auto const target = target_at<Scalar>(index);
        auto const result = probe(target);
        char where[48]{};
        std::snprintf(where, sizeof(where), "walk %+7lld ulp", static_cast<long long>(index - origin));
        if (target > frame.q0) {
            record_walk(result, where, counts, frame.print_points, frame.configuration);
        } else {
            record_behind_start(result, where, counts, frame.print_points, frame.configuration);
        }
    }
}

template <typename Scalar, typename Visit>
void for_each_rung(Scalar q0, Scalar reference, Scalar seed_distance, Visit visit)
{
    // Only commands that still point at the target are recorded. A rung that
    // lands at or behind the start position asks a different question -- the
    // direction test's, not the overshoot test's -- and the two must not be
    // pooled into one population.
    auto const record = [&](Scalar target, char const* rung) {
        if (target > q0) {
            visit(target, rung);
        }
    };
    auto const origin = target_index(reference);
    for (auto exponent = ladder_low_exponent; exponent <= ladder_high_exponent<Scalar>; ++exponent) {
        auto const offset = target_index_t<Scalar>{1} << exponent;
        record(target_at<Scalar>(origin + offset), "ulp");
        record(target_at<Scalar>(origin - offset), "ulp");
    }
    for (std::int32_t decade = 0; decade <= ladder_decades; ++decade) {
        auto const distance = seed_distance * std::pow(Scalar{10}, -static_cast<Scalar>(decade));
        record(reference + distance, "decade");
        record(reference - distance, "decade");
    }
}

template <typename Scalar, typename Probe>
void climb_the_ladder(Scalar q0, Scalar reference, Scalar seed_distance, Probe probe, configuration_counts& counts,
                      bool print_points, char const* configuration)
{
    auto const origin = target_index(reference);
    for_each_rung(q0, reference, seed_distance, [&](Scalar target, char const* rung) {
        auto const result = probe(target);
        ++counts.ladder_swept;
        if (!disagrees(result)) {
            return;
        }
        ++counts.ladder_disagreeing;
        if (print_points) {
            char where[48]{};
            std::snprintf(where, sizeof(where), "ladder %-6s %+lld ulp", rung,
                          static_cast<long long>(target_index(target) - origin));
            print_disagreement(configuration, where, result);
        }
    });
}

}

#endif
