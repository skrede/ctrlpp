#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_BOUNDARY_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_BOUNDARY_H

#include "common.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cstdio>
#include <cstdint>

namespace planner_overshoot_sweep
{

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

inline auto disagrees(verdicts const& result) -> bool
{
    return result.planner_overshoots != result.shape_refuses;
}

inline auto with_shape(verdicts planner, verdicts const& shape) -> verdicts
{
    planner.shape_refuses = shape.shape_refuses;
    planner.shape_error = shape.shape_error;
    planner.shape_error_valid = shape.shape_error_valid;
    return planner;
}

inline void print_disagreement(char const* configuration, char const* where, verdicts const& result)
{
    std::printf("      DISAGREE  %s  %s  planner=%-10s shape=%-7s  reported=%s  shape_error=%s\n", configuration,
                where, result.planner_overshoots ? "overshoots" : "admits",
                result.shape_refuses ? "refuses" : "admits", reason_name(result.reason),
                result.shape_error_valid ? error_name(result.shape_error) : "-");
}

/// Where a verdict turns over, or a statement that it does not.
///
/// The bracket is widened until the verdict holds at its low end and fails at
/// its high end, and both ends are then verified. A verdict that survives the
/// widening in either direction is recorded as such: that is a measurement, not
/// a failure of the search, and it turns out to be the shape of the largest
/// disagreement this sweep finds.
template <typename Scalar>
struct located_boundary
{
    Scalar target{};
    bool located{};
    bool true_across_the_bracket{};
    bool false_across_the_bracket{};
    Scalar bracket_low{};
    Scalar bracket_high{};
};

template <typename Scalar, typename Verdict>
auto widen_bracket(Scalar q0, Scalar seed_distance, Verdict verdict) -> located_boundary<Scalar>
{
    auto low_distance = seed_distance * Scalar{0.5};
    for (std::int32_t attempt = 0; attempt < 60 && !verdict(q0 + low_distance); ++attempt) {
        low_distance *= Scalar{0.5};
    }
    auto high_distance = seed_distance * Scalar{4};
    for (std::int32_t attempt = 0; attempt < 60 && verdict(q0 + high_distance); ++attempt) {
        high_distance *= Scalar{2};
    }

    located_boundary<Scalar> result{};
    result.bracket_low = q0 + low_distance;
    result.bracket_high = q0 + high_distance;
    result.false_across_the_bracket = !verdict(result.bracket_low);
    result.true_across_the_bracket = !result.false_across_the_bracket && verdict(result.bracket_high);
    return result;
}

template <typename Scalar, typename Verdict>
auto bisect_bracket(located_boundary<Scalar> result, Verdict verdict) -> located_boundary<Scalar>
{
    auto low = target_index(result.bracket_low);
    auto high = target_index(result.bracket_high);
    while (high - low > 1) {
        auto const middle = low + (high - low) / 2;
        if (verdict(target_at<Scalar>(middle))) {
            low = middle;
        } else {
            high = middle;
        }
    }
    result.target = target_at<Scalar>(high);
    result.located = true;
    return result;
}

template <typename Scalar, typename Verdict>
auto locate_boundary(Scalar q0, Scalar seed_distance, Verdict verdict) -> located_boundary<Scalar>
{
    auto const bracket = widen_bracket(q0, seed_distance, verdict);
    if (bracket.false_across_the_bracket || bracket.true_across_the_bracket) {
        return bracket;
    }
    return bisect_bracket(bracket, verdict);
}

}

#endif
