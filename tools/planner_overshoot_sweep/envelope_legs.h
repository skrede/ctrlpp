#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ENVELOPE_LEGS_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ENVELOPE_LEGS_H

#include "common.h"
#include "envelope_scan.h"
#include "retarget_windows.h"

#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cstdio>
#include <cstdint>

namespace planner_overshoot_sweep
{

template <typename Scalar>
void record_reason(envelope_scan<Scalar>& scan, Scalar t, bool reversing)
{
    auto const expected = reversing ? online_planner_substitution_reason::reversal_or_overshoot
                                    : online_planner_substitution_reason::none;
    auto const reason = scan.planner.diagnostics().substitution_reason;
    if (reason == expected) {
        return;
    }
    ++scan.counts.unexpected_reason;
    if (scan.print_points) {
        std::printf("      REASON  [%s]  t=%.17g  reported=%s\n", scan.configuration, static_cast<double>(t),
                    reason_name(reason));
    }
}

template <typename Scalar>
void retarget_once(envelope_scan<Scalar>& scan, Scalar t, bool reversing)
{
    auto& planner = scan.planner;
    auto const v_max = scan.limits.v_max;
    planner.reset(Scalar{0});
    auto const driven = planner.update(Scalar{1000} * v_max).has_value();
    auto const q = planner.sample(t).position[0];
    auto const target = reversing ? q - Scalar{1000} * v_max : Scalar{2000} * v_max;
    if (!driven || !planner.update(target).has_value()) {
        std::printf("      REFUSED  [%s]  t=%.17g\n", scan.configuration, static_cast<double>(t));
        return;
    }
    ++scan.counts.retargets;
    record_reason(scan, t, reversing);
    scan_plan(scan, t);
}

template <typename Scalar>
void measure_window(online_planner_3rd<Scalar>& planner, typed_limits<Scalar> const& limits,
                    retarget_window window, bool reversing, bool print_points, envelope_counts& total)
{
    char configuration[128]{};
    std::snprintf(configuration, sizeof(configuration), "v=%.4g a=%.4g j=%.4g %s", static_cast<double>(limits.v_max),
                  static_cast<double>(limits.a_max), static_cast<double>(limits.j_max), window_name(window));
    envelope_counts row{};
    envelope_scan<Scalar> scan{planner, limits, row, configuration, print_points, Scalar{0}};
    auto const dt = static_cast<Scalar>(grid_step);
    auto const span = window_span_of(window, limits, dt);
    for (std::int32_t i = 0; i < instants_per_window; ++i) {
        retarget_once(scan, window_instant(span, dt, i, instants_per_window), reversing);
    }
    std::printf("  [%s]  retargets %lld  reason %lld  |v|>v_max %lld  |a|>a_max %lld  jerk %lld"
                "  allowance %.3e  j_max dt %.3e\n",
                configuration, static_cast<long long>(row.retargets), static_cast<long long>(row.unexpected_reason),
                static_cast<long long>(row.velocity_violations), static_cast<long long>(row.acceleration_violations),
                static_cast<long long>(row.jerk_violations), static_cast<double>(scan.largest_allowance),
                static_cast<double>(limits.j_max * dt));
    accumulate(total, row);
}

inline void report_count(char const* label, std::int64_t count)
{
    std::printf("      %-44s %lld\n", label, static_cast<long long>(count));
}

inline void report_envelope(char const* leg, char const* type, char const* reason_label, envelope_counts const& counts)
{
    std::printf("\n  %s, %s totals\n", leg, type);
    report_count("retargets", counts.retargets);
    report_count(reason_label, counts.unexpected_reason);
    report_count("velocity samples above the limit", counts.velocity_violations);
    report_count("acceleration samples above the limit", counts.acceleration_violations);
    report_count("jerk pairs above the allowance", counts.jerk_violations);
}

template <typename Scalar>
void run_retarget_leg(char const* leg, bool reversing, bool print_points)
{
    std::printf("\n%s, %s\n\n", leg, scalar_name<Scalar>());
    envelope_counts total{};
    auto const windows = reversing ? 2 : 3;
    for (auto const& raw : retarget_limits) {
        auto const limits = typed<Scalar>(raw);
        auto made = online_planner_3rd<Scalar>::create(
            {.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
        if (!made.has_value()) {
            std::printf("  [v=%.4g]  SKIPPED: %s\n", raw.v_max, error_name(made.error()));
            continue;
        }
        for (auto w = 0; w < windows; ++w) {
            measure_window(*made, limits, retarget_windows[w], reversing, print_points, total);
        }
    }
    report_envelope(leg, scalar_name<Scalar>(), reversing ? "reason other than reversal_or_overshoot" : "substituted",
                    total);
}

}

#endif
