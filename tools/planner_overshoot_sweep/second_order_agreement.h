#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SECOND_ORDER_AGREEMENT_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SECOND_ORDER_AGREEMENT_H

#include "common.h"
#include "boundary.h"
#include "band_walk.h"
#include "band_measure.h"

#include "ctrlpp/trajectory/online_planner_2nd.h"

#include <cmath>
#include <cstdio>

namespace planner_overshoot_sweep
{

template <typename Scalar>
struct ramp_state
{
    Scalar q0{};
    Scalar v0{};
    bool reached{};
};

/// Sample the trapezoidal planner part way up its acceleration ramp.
///
/// This planner carries no acceleration state, so every speed on the ramp is a
/// state a caller can present to it, and the reachable speed range is swept
/// rather than probed at one point.
template <typename Scalar>
auto drive_to_ramp(online_planner_2nd<Scalar>& planner, typed_limits<Scalar> const& limits, Scalar q_start,
                   Scalar speed_fraction) -> ramp_state<Scalar>
{
    planner.reset(q_start);

    auto const ramp_time = limits.v_max / limits.a_max;
    auto const drive_time = speed_fraction * ramp_time;
    auto const far_target = q_start + limits.v_max * ramp_time * Scalar{200};

    if (!planner.update(far_target).has_value()) {
        return {};
    }

    auto const point = planner.sample(drive_time);
    ramp_state<Scalar> state{};
    state.q0 = point.position[0];
    state.v0 = point.velocity[0];
    state.reached = (state.v0 > Scalar{0}) && (state.v0 <= limits.v_max) && (state.q0 > Scalar{0});
    return state;
}

template <typename Scalar>
auto report_ramp(ramp_state<Scalar> const& ramp) -> bool
{
    std::printf("      ramp state            q0=%.17g  v0=%.17g\n", static_cast<double>(ramp.q0),
                static_cast<double>(ramp.v0));
    if (!ramp.reached) {
        std::printf("      SKIPPED: no state on the acceleration ramp was reached\n");
    }
    return ramp.reached;
}

template <typename Scalar>
void measure_second_order(online_planner_2nd<Scalar>& planner, typed_limits<Scalar> const& limits,
                          ramp_state<Scalar> const& ramp, char const* configuration, bool print_points,
                          leg_counts& leg)
{
    // The braking distance, mean speed times braking duration, which the planner documents as both its overshoot
    // verdict and its stopping point: a commanded displacement below it has no carry-velocity trapezoid.
    auto const shape_verdict = [&](Scalar target) {
        auto const remaining = target - ramp.q0;
        auto const braking_duration = std::abs(ramp.v0) / limits.a_max;
        auto const braking_distance = std::abs(ramp.v0 * braking_duration / Scalar{2});
        return std::abs(remaining) < braking_distance;
    };
    auto const probe = [&](Scalar target) {
        auto result = planner_verdicts(planner, target);
        result.shape_refuses = shape_verdict(target);
        return result;
    };
    auto const planner_verdict = [&](Scalar target) { return planner_verdicts(planner, target).planner_overshoots; };
    auto const seed_distance = ramp.v0 * ramp.v0 / (Scalar{2} * limits.a_max);
    auto const length_scale = limits.v_max * limits.v_max / (Scalar{2} * limits.a_max);
    configuration_frame<Scalar> const frame{configuration, ramp.q0, seed_distance, length_scale, print_points};
    measure_configuration(frame, planner_verdict, shape_verdict, probe, leg);
}

template <typename Scalar>
void sweep_second_order_configuration(limit_set const& raw, double q_start, double speed_fraction,
                                      bool print_points, leg_counts& leg)
{
    ++leg.configurations;
    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g q0~%.4g v0/v=%.2f", raw.v_max, raw.a_max,
                  q_start, speed_fraction);
    auto const limits = typed<Scalar>(raw);
    auto made = online_planner_2nd<Scalar>::create({.v_max = limits.v_max, .a_max = limits.a_max});
    std::printf("  [%s]\n", configuration);
    if (!made.has_value()) {
        std::printf("      SKIPPED: the limits were rejected (%s)\n", error_name(made.error()));
        return;
    }
    auto const ramp =
        drive_to_ramp(*made, limits, static_cast<Scalar>(q_start), static_cast<Scalar>(speed_fraction));
    if (!report_ramp(ramp)) {
        return;
    }
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g q0=%.17g v0=%.17g", raw.v_max, raw.a_max,
                  static_cast<double>(ramp.q0), static_cast<double>(ramp.v0));
    measure_second_order(*made, limits, ramp, configuration, print_points, leg);
}

}

#endif
