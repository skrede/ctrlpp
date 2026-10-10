#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_AGREEMENT_LEGS_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_AGREEMENT_LEGS_H

#include "common.h"
#include "band_measure.h"
#include "third_order_agreement.h"
#include "second_order_agreement.h"

#include <cstdio>
#include <limits>

namespace planner_overshoot_sweep
{

/// Limit sets spanning several decades in every limit, chosen so that both
/// branches of the jerk-limited stopping-distance computation are covered: the
/// deceleration reaches the acceleration limit when the cruise speed exceeds
/// the squared acceleration limit over the jerk limit, and does not otherwise.
/// Which branch each set takes is printed with its result rather than asserted
/// here.
inline constexpr limit_set third_order_limits[] = {
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
inline constexpr double axis_positions[] = {0.0, 1e2, 1e4};

inline constexpr limit_set second_order_limits[] = {
    {5.0, 10.0, 0.0}, {1.0, 5.0, 0.0}, {0.25, 1.0, 0.0}, {40.0, 20.0, 0.0}, {1e3, 1e4, 0.0},
};

inline constexpr double second_order_speed_fractions[] = {0.1, 0.35, 0.7, 0.95};

inline void report_leg(char const* name, leg_counts const& counts)
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
    std::printf("      walk targets at or behind the start, skipped                 %lld\n", counts.walk_behind_start);
    std::printf("      behind-start targets not reported as reversal_or_overshoot   %lld\n",
                counts.behind_start_not_reversal);
}

template <typename Scalar>
void run_third_order_leg(bool print_points)
{
    std::printf("Third-order planner: its overshoot test against the double-S constructor's own\n");
    std::printf("displacement rejection, on exactly the command the planner hands it.\n\n");
    leg_counts third{};
    for (auto const& limits : third_order_limits) {
        for (auto const position : axis_positions) {
            sweep_third_order_configuration<Scalar>(limits, position, print_points, third);
        }
    }
    report_leg("Third-order", third);
}

template <typename Scalar>
void run_second_order_leg(bool print_points)
{
    std::printf("\n\nSecond-order planner: its overshoot test against the braking distance\n");
    std::printf("formula it documents for its decision and its stopping point. The double-S\n");
    std::printf("constructor is not reachable from this planner, so no rejection of its own\n");
    std::printf("exists to compare against.\n\n");
    leg_counts second{};
    for (auto const& limits : second_order_limits) {
        for (auto const position : axis_positions) {
            for (auto const fraction : second_order_speed_fractions) {
                sweep_second_order_configuration<Scalar>(limits, position, fraction, print_points, second);
            }
        }
    }
    report_leg("Second-order", second);
}

template <typename Scalar>
void run_agreement_legs(bool print_points)
{
    std::printf("%s, epsilon = %.6e\n\n", scalar_name<Scalar>(),
                static_cast<double>(std::numeric_limits<Scalar>::epsilon()));
    run_third_order_leg<Scalar>(print_points);
    run_second_order_leg<Scalar>(print_points);
}

}

#endif
