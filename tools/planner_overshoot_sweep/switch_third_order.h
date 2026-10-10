#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_THIRD_ORDER_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_THIRD_ORDER_H

#include "common.h"
#include "switch_scan.h"
#include "switch_planners.h"
#include "switch_population.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <type_traits>

namespace planner_overshoot_sweep
{

inline constexpr double planner_resolution_multiples[] = {0.25, 0.5, 1.0, 2.0, 4.0};

template <typename Scalar>
struct third_order_command
{
    typed_limits<Scalar> limits;
    Scalar target;
    Scalar cruise_fraction;
    Scalar ramp_fraction;
};

template <typename Scalar>
void scan_third_order(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits, Scalar target,
                      std::vector<Scalar> const& centers, Scalar t_r, bool closed_form, planner_leg& leg)
{
    char configuration[200];
    std::snprintf(configuration, sizeof(configuration), "v=%a a=%a j=%a h=%a t_r=%a", static_cast<double>(limits.v_max),
                  static_cast<double>(limits.a_max), static_cast<double>(limits.j_max), static_cast<double>(target),
                  static_cast<double>(t_r));
    switch_scan<Scalar> scan{leg.counts, configuration, limits.v_max, limits.a_max, switch_evidence::acceleration_changes,
                             closed_form, leg.print_points};
    scan_switches(scan, centers, planner_domain<Scalar>(), [&](Scalar t) { return planner.sample(t); });
}

template <typename Scalar>
auto plan_unsubstituted(online_planner_3rd<Scalar>& planner, Scalar target) -> bool
{
    return planner.update(target).has_value()
           && planner.diagnostics().substitution_reason == online_planner_substitution_reason::none;
}

template <typename Scalar>
auto drive_and_sample(online_planner_3rd<Scalar>& planner, Scalar target, Scalar t) -> trajectory_point<Scalar, 1>
{
    planner.reset(Scalar{0});
    (void)planner.update(target);
    return planner.sample(t);
}

template <typename Scalar>
struct nulling_end
{
    Scalar duration;
    Scalar q;
    Scalar v;
};

// @cite biagiotti2009 -- Sec. 4.6.1: the acceleration is first brought to zero at the jerk limit, then a double-S
// carries the speed it leaves. The planner's own trim of that phase is private and not mirrored here.
template <typename Scalar>
auto null_acceleration(typed_limits<Scalar> const& limits, trajectory_point<Scalar, 1> const& state) -> nulling_end<Scalar>
{
    auto const q = state.position[0];
    auto const v = state.velocity[0];
    auto const a = state.acceleration[0];
    auto const T = std::abs(a) / limits.j_max;
    auto const j = a > Scalar{0} ? -limits.j_max : limits.j_max;
    auto const v_end = v + a * T + Scalar{0.5} * j * T * T;
    auto const q_end = q + v * T + Scalar{0.5} * a * T * T + j * T * T * T / Scalar{6};
    return {T, q_end, std::copysign(std::min(std::abs(v_end), limits.v_max), v_end)};
}

template <typename Scalar>
void measure_from_cruise(online_planner_3rd<Scalar>& planner, third_order_command<Scalar> const& c,
                         double_s_trajectory<Scalar> const& first, planner_leg& leg)
{
    ++leg.counts.configurations;
    auto const p = first.phase_durations();
    auto const t_r = p[0] + p[1] + p[2] + c.cruise_fraction * p[3];
    auto const state = drive_and_sample(planner, c.target, t_r);
    auto const second = carry_shape(c.limits, state.position[0], Scalar{2} * c.target, state.velocity[0]);
    if (state.acceleration[0] != Scalar{0} || !second.has_value() || !plan_unsubstituted(planner, Scalar{2} * c.target)) {
        ++leg.counts.not_measured;
        return;
    }
    scan_third_order(planner, c.limits, Scalar{2} * c.target, planner_switches(*second, t_r), t_r, false, leg);
}

template <typename Scalar>
void measure_from_ramp(online_planner_3rd<Scalar>& planner, third_order_command<Scalar> const& c,
                       double_s_trajectory<Scalar> const& first, planner_leg& leg)
{
    ++leg.counts.configurations;
    auto const p = first.phase_durations();
    auto const t_r = c.ramp_fraction * (p[0] + p[1] + p[2]);
    auto const state = drive_and_sample(planner, c.target, t_r);
    auto const end = null_acceleration(c.limits, state);
    auto const second = carry_shape(c.limits, end.q, Scalar{2} * c.target, end.v);
    if (!second.has_value() || !plan_unsubstituted(planner, Scalar{2} * c.target)) {
        ++leg.counts.not_measured;
        return;
    }
    auto centers = planner_switches(*second, t_r + end.duration);
    if (end.duration > Scalar{0}) {
        centers.push_back(t_r + end.duration);
    }
    scan_third_order(planner, c.limits, Scalar{2} * c.target, centers, t_r, true, leg);
}

template <typename Scalar>
void measure_third_order(online_planner_3rd<Scalar>& planner, third_order_command<Scalar> const& c, planner_leg& leg)
{
    ++leg.counts.configurations;
    planner.reset(Scalar{0});
    auto const first = carry_shape(c.limits, Scalar{0}, c.target, Scalar{0});
    if (!first.has_value() || !plan_unsubstituted(planner, c.target)) {
        leg.counts.configurations += 2;
        leg.counts.not_measured += 3;
        return;
    }
    scan_third_order(planner, c.limits, c.target, planner_switches(*first, Scalar{0}), Scalar{0}, false, leg);
    measure_from_cruise(planner, c, *first, leg);
    measure_from_ramp(planner, c, *first, leg);
}

template <typename Scalar>
void measure_third_order_limits(typed_limits<Scalar> const& limits, fraction_source<Scalar>& source, planner_leg& leg)
{
    auto made = online_planner_3rd<Scalar>::create({.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
    for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
        auto const target = sigma * Scalar{4} * limits.v_max * acceleration_stretch(limits);
        auto const cruise_fraction = source.unit_fraction();
        auto const command = third_order_command<Scalar>{limits, target, cruise_fraction, source.unit_fraction()};
        if (!made.has_value()) {
            leg.counts.configurations += 3;
            leg.counts.not_measured += 3;
            continue;
        }
        measure_third_order(*made, command, leg);
    }
}

template <typename Scalar>
void run_third_order_switch_leg(bool print_points)
{
    auto const seed = switch_seed<Scalar>(switch_leg::online_planner_3rd);
    fraction_source<Scalar> source{seed};
    planner_leg leg{{}, print_points};
    for (auto const& raw : retarget_limits) {
        for (auto const c : planner_resolution_multiples) {
            auto limits = typed<Scalar>(raw);
            limits.j_max = limits.a_max / (static_cast<Scalar>(c) * resolution_at(acceleration_stretch(limits)));
            measure_third_order_limits(limits, source, leg);
        }
    }
    // Limits on which an earlier random search sampled a speed above v_max inside the constant-acceleration phase.
    if constexpr (std::is_same_v<Scalar, float>) {
        measure_third_order_limits(typed_limits<float>{18.7756824f, 0.0103791943f, 136.34903f}, source, leg);
    }
    report_switch_block("online_planner_3rd", scalar_name<Scalar>(), seed, leg.counts);
}

}

#endif
