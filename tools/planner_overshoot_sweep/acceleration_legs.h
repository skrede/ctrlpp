#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ACCELERATION_LEGS_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ACCELERATION_LEGS_H

#include "common.h"
#include "boundary.h"
#include "envelope_legs.h"
#include "retarget_windows.h"
#include "acceleration_probe.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iterator>

namespace planner_overshoot_sweep
{

// A seed for the bracket only: the cruise leg's stopping distance from the speed the nulling phase ends at.
template <typename Scalar>
auto acceleration_seed(typed_limits<Scalar> const& limits, trajectory_point<Scalar, 1> const& point) -> Scalar
{
    auto const a0 = point.acceleration[0];
    auto const v_peak = point.velocity[0] + a0 * std::abs(a0) / (Scalar{2} * limits.j_max);
    return v_peak * v_peak / (Scalar{2} * limits.a_max) + limits.a_max * v_peak / (Scalar{2} * limits.j_max);
}

template <typename Scalar>
void describe_state(char (&out)[192], typed_limits<Scalar> const& limits, retarget_window window, Scalar t,
                    trajectory_point<Scalar, 1> const& point)
{
    std::snprintf(out, sizeof(out), "v=%.4g a=%.4g j=%.4g %s t=%.9g q0=%.17g v0=%.17g a0=%.17g",
                  static_cast<double>(limits.v_max), static_cast<double>(limits.a_max),
                  static_cast<double>(limits.j_max), window_name(window), static_cast<double>(t),
                  static_cast<double>(point.position[0]), static_cast<double>(point.velocity[0]),
                  static_cast<double>(point.acceleration[0]));
}

template <typename Scalar>
void report_unlocated(char const* configuration, located_boundary<Scalar> const& boundary, Scalar q0)
{
    std::printf("  [%s]  BOUNDARY NOT LOCATED: the overshoot verdict %s across remaining distances from %.6e to %.6e\n",
                configuration, boundary.true_across_the_bracket ? "holds" : "fails",
                static_cast<double>(boundary.bracket_low - q0), static_cast<double>(boundary.bracket_high - q0));
}

inline void report_acceleration_state(char const* configuration, double boundary, acceleration_counts const& row)
{
    std::printf("  [%s]  boundary %.17g  probes %lld  shape refusals %lld  substitutions at or above %lld\n",
                configuration, boundary, static_cast<long long>(row.probes),
                static_cast<long long>(row.shape_refusals), static_cast<long long>(row.substitutions_above));
}

template <typename Scalar>
void measure_acceleration_state(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits,
                                retarget_window window, Scalar t, bool print_points, acceleration_counts& counts)
{
    ++counts.states;
    auto const point = drive_to_instant(planner, limits.v_max, t);
    char configuration[192]{};
    describe_state(configuration, limits, window, t, point);
    auto const overshoots = [&](Scalar target) {
        return reason_after_retarget(planner, limits.v_max, t, target)
               == online_planner_substitution_reason::reversal_or_overshoot;
    };
    auto const seed_distance = acceleration_seed(limits, point);
    auto const boundary = locate_boundary(point.position[0], seed_distance, overshoots);
    if (!boundary.located) {
        report_unlocated(configuration, boundary, point.position[0]);
        return;
    }
    ++counts.boundaries_located;
    acceleration_counts row{};
    acceleration_probe<Scalar> probe{planner, limits.v_max, t, boundary.target, configuration, print_points, row};
    probe_around_boundary(probe, point.position[0], seed_distance);
    report_acceleration_state(configuration, static_cast<double>(boundary.target), row);
    accumulate(counts, row);
}

template <typename Scalar>
void measure_acceleration_windows(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits,
                                  bool print_points, acceleration_counts& counts)
{
    auto const dt = static_cast<Scalar>(grid_step);
    for (auto const window : acceleration_windows) {
        auto const span = window_span_of(window, limits, dt);
        for (std::int32_t i = 0; i < instants_per_acceleration_window; ++i) {
            auto const t = window_instant(span, dt, i, instants_per_acceleration_window);
            measure_acceleration_state(planner, limits, window, t, print_points, counts);
        }
    }
}

inline void report_acceleration(char const* type, acceleration_counts const& counts)
{
    std::printf("\n  Third-order from nonzero acceleration, %s totals\n", type);
    report_count("retarget states", counts.states);
    report_count("boundaries located", counts.boundaries_located);
    report_count("shape refusals reported", counts.shape_refusals);
    report_count("substitutions at or above the boundary", counts.substitutions_above);
}

template <typename Scalar>
void run_acceleration_leg(bool print_points)
{
    std::printf("\nThird-order from nonzero acceleration, %s\n\n", scalar_name<Scalar>());
    acceleration_counts total{};
    for (auto const& raw : retarget_limits) {
        auto const limits = typed<Scalar>(raw);
        auto made = online_planner_3rd<Scalar>::create(
            {.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
        if (!made.has_value()) {
            std::printf("  [v=%.4g]  SKIPPED: %s\n", raw.v_max, error_name(made.error()));
            total.states += std::ssize(acceleration_windows) * instants_per_acceleration_window;
            continue;
        }
        measure_acceleration_windows(*made, limits, print_points, total);
    }
    report_acceleration(scalar_name<Scalar>(), total);
}

inline void run_acceleration_legs(bool print_points)
{
    std::printf("\n\nThird-order planner retargeted forward from states in its acceleration-ramp and\n");
    std::printf("acceleration-nulling windows: each state's overshoot boundary located on the planner's own\n");
    std::printf("report, then walked and laddered for carry-velocity shape refusals and for substitutions\n");
    std::printf("at or above the boundary.\n");
    run_acceleration_leg<double>(print_points);
    run_acceleration_leg<float>(print_points);
}

}

#endif
