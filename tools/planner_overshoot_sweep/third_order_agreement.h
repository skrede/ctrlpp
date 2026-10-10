#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_THIRD_ORDER_AGREEMENT_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_THIRD_ORDER_AGREEMENT_H

#include "common.h"
#include "boundary.h"
#include "band_walk.h"
#include "band_measure.h"

#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cstdio>

namespace planner_overshoot_sweep
{

/// The state the planner reports part way through its own cruise phase.
///
/// Reported by the planner rather than assembled here. The acceleration is
/// required to be exactly zero, which is what keeps the acceleration-nulling
/// phase out of the command that follows: with it out, the command reaches the
/// overshoot test carrying the position and speed printed here, and the shape
/// constructor can be handed the same two numbers.
template <typename Scalar>
struct cruise_state
{
    Scalar q0{};
    Scalar v0{};
    Scalar a0{};
    bool reached{};
};

template <typename Scalar>
auto drive_to_cruise(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits, Scalar q_start)
    -> cruise_state<Scalar>
{
    planner.reset(q_start);

    auto const acceleration_stretch = limits.a_max / limits.j_max + limits.v_max / limits.a_max;
    auto const drive_time = Scalar{4} * acceleration_stretch;
    auto const far_target = q_start + limits.v_max * drive_time * Scalar{20};

    if (!planner.update(far_target).has_value()) {
        return {};
    }

    auto const point = planner.sample(drive_time);
    cruise_state<Scalar> state{};
    state.q0 = point.position[0];
    state.v0 = point.velocity[0];
    state.a0 = point.acceleration[0];
    state.reached = (state.a0 == Scalar{0}) && (state.v0 > Scalar{0}) && (state.q0 > Scalar{0});
    return state;
}

template <typename Scalar>
auto report_cruise(cruise_state<Scalar> const& cruise, typed_limits<Scalar> const& limits) -> bool
{
    std::printf("      cruise state          q0=%.17g  v0=%.17g  a0=%.17g\n", static_cast<double>(cruise.q0),
                static_cast<double>(cruise.v0), static_cast<double>(cruise.a0));
    if (!cruise.reached) {
        std::printf("      SKIPPED: no cruise state with an exactly zero acceleration was"
                    " reached\n");
        return false;
    }
    std::printf("      cruise speed against the velocity limit   %+lld ulp\n",
                static_cast<long long>(ulps_between(limits.v_max, cruise.v0)));
    std::printf("      deceleration reaches the acceleration limit   %s\n",
                (cruise.v0 > limits.a_max * limits.a_max / limits.j_max) ? "yes" : "no");
    return true;
}

/// The constructor's verdict, called on exactly the configuration the planner
/// hands it when it builds a carry-velocity profile: the same start position,
/// the same carried velocity, the same limits and a terminal velocity of zero.
template <typename Scalar>
auto double_s_verdicts(typed_limits<Scalar> const& limits, cruise_state<Scalar> const& cruise, Scalar target)
    -> verdicts
{
    auto const built = double_s_trajectory<Scalar>::create({.q0 = cruise.q0,
                                                            .q1 = target,
                                                            .v_max = limits.v_max,
                                                            .a_max = limits.a_max,
                                                            .j_max = limits.j_max,
                                                            .v0 = cruise.v0,
                                                            .v1 = Scalar{0}});
    verdicts result{};
    if (built.has_value()) {
        return result;
    }
    result.shape_refuses = true;
    result.shape_error = built.error();
    result.shape_error_valid = true;
    return result;
}

template <typename Scalar>
void measure_third_order(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits,
                         cruise_state<Scalar> const& cruise, char const* configuration, bool print_points,
                         leg_counts& leg)
{
    auto const shape_probe = [&](Scalar target) { return double_s_verdicts(limits, cruise, target); };
    auto const probe = [&](Scalar target) { return with_shape(planner_verdicts(planner, target), shape_probe(target)); };
    auto const planner_verdict = [&](Scalar target) { return planner_verdicts(planner, target).planner_overshoots; };
    auto const shape_verdict = [&](Scalar target) { return shape_probe(target).shape_refuses; };
    // A seed for the bracket only. It decides nothing: every bracket end is
    // widened until the verdict holds there and is then verified.
    auto const seed_distance = cruise.v0 * cruise.v0 / (Scalar{2} * limits.a_max)
                               + limits.a_max * cruise.v0 / (Scalar{2} * limits.j_max);
    // The direction test's length floor is a counted multiple of epsilon times the stopping distance from full
    // speed; the remaining distance in that scale shows the floor clause did not select a substitution here.
    auto const length_scale = limits.v_max * limits.v_max / (Scalar{2} * limits.a_max);
    configuration_frame<Scalar> const frame{configuration, cruise.q0, seed_distance, length_scale, print_points};
    measure_configuration(frame, planner_verdict, shape_verdict, probe, leg);
}

template <typename Scalar>
void sweep_third_order_configuration(limit_set const& raw, double q_start, bool print_points, leg_counts& leg)
{
    ++leg.configurations;
    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g j=%.4g q0~%.4g", raw.v_max, raw.a_max,
                  raw.j_max, q_start);
    auto const limits = typed<Scalar>(raw);
    auto made = online_planner_3rd<Scalar>::create({.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
    std::printf("  [%s]\n", configuration);
    if (!made.has_value()) {
        std::printf("      SKIPPED: the limits were rejected (%s)\n", error_name(made.error()));
        return;
    }
    auto const cruise = drive_to_cruise(*made, limits, static_cast<Scalar>(q_start));
    if (!report_cruise(cruise, limits)) {
        return;
    }
    // The per-point lines carry the speed the command was planned from as well
    // as the limits, so a disagreeing point is fully identified by its own line.
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g j=%.4g q0=%.17g v0=%.17g", raw.v_max,
                  raw.a_max, raw.j_max, static_cast<double>(cruise.q0), static_cast<double>(cruise.v0));
    measure_third_order(*made, limits, cruise, configuration, print_points, leg);
}

}

#endif
