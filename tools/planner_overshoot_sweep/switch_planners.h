#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_PLANNERS_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_PLANNERS_H

#include "common.h"
#include "switch_scan.h"
#include "switch_population.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
#include <cstdint>

namespace planner_overshoot_sweep
{

struct planner_leg
{
    switch_counts counts;
    bool print_points;
};

template <typename Scalar>
auto planner_domain() -> instant_span<Scalar>
{
    return {Scalar{0}, std::numeric_limits<Scalar>::max()};
}

template <typename Scalar>
auto carry_shape(typed_limits<Scalar> const& limits, Scalar q0, Scalar q1, Scalar v0)
    -> ctrlpp::expected<double_s_trajectory<Scalar>, trajectory_error>
{
    return double_s_trajectory<Scalar>::create(
        {.q0 = q0, .q1 = q1, .v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max, .v0 = v0, .v1 = Scalar{0}});
}

// The planner appends the shape's nonzero phases and switches at their running sum from the plan's start.
template <typename Scalar>
auto planner_switches(double_s_trajectory<Scalar> const& shape, Scalar start) -> std::vector<Scalar>
{
    std::vector<Scalar> centers{start + shape.duration()};
    auto cumulative = Scalar{0};
    for (auto const phase : shape.phase_durations()) {
        cumulative += phase;
        if (phase > Scalar{0}) {
            centers.push_back(start + cumulative);
        }
    }
    return centers;
}

// The second-order leg also runs a pair with full-width significands: on the short ones above, v_max - v_r is
// often exact, which hides the rounding this leg looks for.
inline constexpr limit_set full_mantissa_limits{0.006846007648171031, 0.20062093968092909, 0.0};

template <typename Scalar>
struct same_direction_retarget
{
    Scalar h1;
    Scalar fraction;
    Scalar h2;
};

template <typename Scalar>
auto draw_retarget(fraction_source<Scalar>& source, typed_limits<Scalar> const& limits, Scalar sigma) -> same_direction_retarget<Scalar>
{
    auto const scale = limits.v_max * limits.v_max / limits.a_max;
    auto const h1 = sigma * source.magnitude() * scale;
    auto const fraction = source.unit_fraction();
    return {h1, fraction, h1 + sigma * source.magnitude() * scale};
}

// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14): the ramp reaches v_max exactly when the triangular peak does.
template <typename Scalar>
auto reaches_limit(typed_limits<Scalar> const& limits, Scalar remaining, Scalar v_r) -> bool
{
    return std::sqrt((Scalar{2} * limits.a_max * std::abs(remaining) + v_r * v_r) / Scalar{2}) >= limits.v_max;
}

// The retarget is taken on the first half of the first plan's ramp, where v_r < v_max / 2; the new ramp then
// reaches v_max after (v_max - |v_r|) / a_max.
template <typename Scalar>
void measure_second_order(online_planner_2nd<Scalar>& planner, typed_limits<Scalar> const& limits,
                          same_direction_retarget<Scalar> const& r, planner_leg& leg)
{
    ++leg.counts.configurations;
    planner.reset(Scalar{0});
    auto const first = planner.update(r.h1).has_value();
    auto const t_r = r.fraction * limits.v_max / limits.a_max / Scalar{2};
    auto const state = planner.sample(t_r);
    auto const v_r = state.velocity[0];
    if (!first || !planner.update(r.h2).has_value() || planner.diagnostics().substitution_reason != online_planner_substitution_reason::none
        || !(std::abs(v_r) < limits.v_max) || !reaches_limit(limits, r.h2 - state.position[0], v_r)) {
        ++leg.counts.not_measured;
        return;
    }
    char configuration[200];
    std::snprintf(configuration, sizeof(configuration), "v=%a a=%a h1=%a t_r=%a h2=%a", static_cast<double>(limits.v_max),
                  static_cast<double>(limits.a_max), static_cast<double>(r.h1), static_cast<double>(t_r), static_cast<double>(r.h2));
    switch_scan<Scalar> scan{leg.counts,  configuration, limits.v_max, limits.a_max, switch_evidence::limit_then_zero,
                             false,       leg.print_points};
    auto const t_switch = t_r + (limits.v_max - std::abs(v_r)) / limits.a_max;
    scan_switches(scan, {t_switch}, planner_domain<Scalar>(), [&](Scalar t) { return planner.sample(t); });
}

template <typename Scalar>
void measure_second_order_limits(limit_set const& raw, fraction_source<Scalar>& source, planner_leg& leg)
{
    auto const limits = typed<Scalar>(raw);
    auto made = online_planner_2nd<Scalar>::create({.v_max = limits.v_max, .a_max = limits.a_max});
    for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
        for (std::int32_t k = 0; k < seeded_pairs; ++k) {
            auto const retarget = draw_retarget(source, limits, sigma);
            if (made.has_value()) {
                measure_second_order(*made, limits, retarget, leg);
            }
        }
    }
}

template <typename Scalar>
void run_second_order_switch_leg(bool print_points)
{
    auto const seed = switch_seed<Scalar>(switch_leg::online_planner_2nd);
    fraction_source<Scalar> source{seed};
    planner_leg leg{{}, print_points};
    for (auto const& raw : retarget_limits) {
        measure_second_order_limits(raw, source, leg);
    }
    measure_second_order_limits(full_mantissa_limits, source, leg);
    report_switch_block("online_planner_2nd", scalar_name<Scalar>(), seed, leg.counts);
}

}

#endif
