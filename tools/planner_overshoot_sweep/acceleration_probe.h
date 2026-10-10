#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ACCELERATION_PROBE_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ACCELERATION_PROBE_H

#include "common.h"
#include "band_walk.h"
#include "retarget_windows.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cstdio>
#include <cstdint>

namespace planner_overshoot_sweep
{

inline constexpr std::int32_t instants_per_acceleration_window = 5;

inline constexpr retarget_window acceleration_windows[] = {retarget_window::acceleration_ramp,
                                                           retarget_window::acceleration_nulling};

struct acceleration_counts
{
    std::int64_t states;
    std::int64_t boundaries_located;
    std::int64_t probes;
    std::int64_t shape_refusals;
    std::int64_t substitutions_above;
};

inline void accumulate(acceleration_counts& total, acceleration_counts const& row)
{
    total.states += row.states;
    total.boundaries_located += row.boundaries_located;
    total.probes += row.probes;
    total.shape_refusals += row.shape_refusals;
    total.substitutions_above += row.substitutions_above;
}

// Every probe drives the planner afresh from rest to the same instant, so each command meets the identical state.
template <typename Scalar>
auto drive_to_instant(online_planner_3rd<Scalar>& planner, Scalar v_max, Scalar t) -> trajectory_point<Scalar, 1>
{
    planner.reset(Scalar{0});
    if (!planner.update(Scalar{1000} * v_max).has_value()) {
        std::printf("      REFUSED  drive  v_max=%.17g\n", static_cast<double>(v_max));
    }
    return planner.sample(t);
}

template <typename Scalar>
auto reason_after_retarget(online_planner_3rd<Scalar>& planner, Scalar v_max, Scalar t, Scalar target)
    -> online_planner_substitution_reason
{
    drive_to_instant(planner, v_max, t);
    if (!planner.update(target).has_value()) {
        std::printf("      REFUSED  retarget  t=%.17g  target=%.17g\n", static_cast<double>(t),
                    static_cast<double>(target));
    }
    return planner.diagnostics().substitution_reason;
}

template <typename Scalar>
struct acceleration_probe
{
    online_planner_3rd<Scalar>& planner;
    Scalar v_max;
    Scalar t;
    Scalar boundary;
    char const* configuration;
    bool print_points;
    acceleration_counts& counts;
};

// A refusal by the carry-velocity shape is a command the planner's own displacement test admitted, and a substitution
// at or above the planner's boundary is a brake to rest with room to stop.
template <typename Scalar>
void probe_target(acceleration_probe<Scalar>& probe, Scalar target, char const* where)
{
    auto const reason = reason_after_retarget(probe.planner, probe.v_max, probe.t, target);
    auto const refused = reason == online_planner_substitution_reason::carry_velocity_shape_unavailable;
    auto const substituted = target >= probe.boundary && reason != online_planner_substitution_reason::none;
    ++probe.counts.probes;
    probe.counts.shape_refusals += refused ? 1 : 0;
    probe.counts.substitutions_above += substituted ? 1 : 0;
    if (probe.print_points && (refused || substituted)) {
        std::printf("      DISAGREE  [%s]  %s %+lld ulp  reported=%s\n", probe.configuration, where,
                    static_cast<long long>(target_index(target) - target_index(probe.boundary)), reason_name(reason));
    }
}

template <typename Scalar>
void probe_around_boundary(acceleration_probe<Scalar>& probe, Scalar q0, Scalar seed_distance)
{
    auto const origin = target_index(probe.boundary);
    for (auto index = origin - walk_margin_ulps; index <= origin + walk_margin_ulps; ++index) {
        probe_target(probe, target_at<Scalar>(index), "walk");
    }
    for_each_rung(q0, probe.boundary, seed_distance,
                  [&](Scalar target, char const* rung) { probe_target(probe, target, rung); });
}

}

#endif
