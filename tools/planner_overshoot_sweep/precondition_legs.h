#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_PRECONDITION_LEGS_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_PRECONDITION_LEGS_H

#include "common.h"
#include "envelope_legs.h"
#include "retarget_windows.h"
#include "precondition_build.h"

#include "ctrlpp/trajectory/trajectory_types.h"

#include <cstdio>
#include <cstddef>
#include <cstdint>

namespace planner_overshoot_sweep
{

template <typename Shape, typename Scalar>
void measure_command(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command,
                     char const* configuration, bool print_points, precondition_counts& counts)
{
    ++counts.commands;
    char outcome[96]{};
    auto const built = build_shape<Shape>(limits, command);
    if (built.has_value()) {
        auto const peak = record_build(limits, command, *built, counts);
        std::snprintf(outcome, sizeof(outcome), "built  peak %.17g", static_cast<double>(peak));
    } else {
        record_refusal(limits, command, built.error(), counts);
        std::snprintf(outcome, sizeof(outcome), "refused  %s", error_name(built.error()));
    }
    if (print_points) {
        print_command(configuration, command, outcome);
    }
}

inline void print_precondition_row(char const* configuration, precondition_counts const& row)
{
    std::printf("  [%s]  commands %lld  builds %lld", configuration, static_cast<long long>(row.commands),
                static_cast<long long>(row.builds));
    for (std::size_t e = 0; e < error_count; ++e) {
        if (row.refusals[e] != 0) {
            std::printf("  %s %lld", error_name(static_cast<trajectory_error>(e)),
                        static_cast<long long>(row.refusals[e]));
        }
    }
    std::printf("  largest peak / v_max - 1 %+.3e\n", row.largest_peak_ratio - 1.0);
}

template <typename Shape, typename Scalar>
void measure_limit_set(typed_limits<Scalar> const& limits, char const* shape_name, bool print_points,
                       precondition_counts& total)
{
    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "%s %s v=%.4g a=%.4g j=%.4g", shape_name,
                  scalar_name<Scalar>(), static_cast<double>(limits.v_max), static_cast<double>(limits.a_max),
                  static_cast<double>(limits.j_max));
    precondition_counts row{};
    for (auto step = lowest_speed_step; step <= highest_speed_step; ++step) {
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            for (auto const final_role : {false, true}) {
                for (auto const knob : knob_values) {
                    boundary_command<Scalar> const command{step, stepped_speed(limits.v_max, step), sigma,
                                                           final_role, knob};
                    measure_command<Shape>(limits, command, configuration, print_points, row);
                }
            }
        }
    }
    print_precondition_row(configuration, row);
    accumulate(total, row);
}

inline void report_precondition(char const* shape_name, char const* type, precondition_counts const& counts)
{
    std::printf("\n  Boundary-velocity precondition, %s, %s totals\n", shape_name, type);
    report_count("commands", counts.commands);
    report_count("built", counts.builds);
    for (std::size_t e = 0; e < error_count; ++e) {
        if (counts.refusals[e] != 0) {
            std::printf("      refused, %-35s %lld\n", error_name(static_cast<trajectory_error>(e)),
                        static_cast<long long>(counts.refusals[e]));
        }
    }
    report_count("builds with a peak above the limit", counts.peak_above_limit);
    report_count("refusals at or below the limit", counts.refused_at_or_below);
    report_count("refused within the allowance above the limit", counts.refused_within_allowance);
    report_count("admissions beyond the allowance", counts.admitted_beyond_allowance);
    report_count("clamped ends reported exactly at the limit", counts.clamped_exact);
    report_count("clamped ends reported off the limit", counts.clamped_off);
}

template <typename Shape>
void run_precondition_leg(char const* shape_name, bool print_points)
{
    using Scalar = typename Shape::scalar_type;
    std::printf("\nBoundary-velocity precondition, %s, %s\n\n", shape_name, scalar_name<Scalar>());
    precondition_counts total{};
    for (auto const& raw : retarget_limits) {
        measure_limit_set<Shape>(typed<Scalar>(raw), shape_name, print_points, total);
    }
    report_precondition(shape_name, scalar_name<Scalar>(), total);
}

}

#endif
