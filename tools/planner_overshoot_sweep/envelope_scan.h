#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ENVELOPE_SCAN_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_ENVELOPE_SCAN_H

#include "common.h"
#include "retarget_windows.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdint>
#include <algorithm>

namespace planner_overshoot_sweep
{

struct envelope_counts
{
    std::int64_t retargets;
    std::int64_t unexpected_reason;
    std::int64_t velocity_violations;
    std::int64_t acceleration_violations;
    std::int64_t jerk_violations;
};

inline void accumulate(envelope_counts& total, envelope_counts const& row)
{
    total.retargets += row.retargets;
    total.unexpected_reason += row.unexpected_reason;
    total.velocity_violations += row.velocity_violations;
    total.acceleration_violations += row.acceleration_violations;
    total.jerk_violations += row.jerk_violations;
}

template <typename Scalar>
struct envelope_scan
{
    online_planner_3rd<Scalar>& planner;
    typed_limits<Scalar> const& limits;
    envelope_counts& counts;
    char const* configuration;
    bool print_points;
    Scalar largest_allowance;
};

template <typename Scalar>
void report_violation(envelope_scan<Scalar> const& scan, char const* what, Scalar t, Scalar value)
{
    if (scan.print_points) {
        std::printf("      VIOLATION  [%s]  %s  t=%.17g  value=%.17g\n", scan.configuration, what,
                    static_cast<double>(t), static_cast<double>(value));
    }
}

template <typename Scalar>
void record_sample(envelope_scan<Scalar>& scan, Scalar t, trajectory_point<Scalar, 1> const& point)
{
    auto const v = std::abs(point.velocity[0]);
    auto const a = std::abs(point.acceleration[0]);
    if (v > scan.limits.v_max) {
        ++scan.counts.velocity_violations;
        report_violation(scan, "|v| above v_max", t, v);
    }
    if (a > scan.limits.a_max) {
        ++scan.counts.acceleration_violations;
        report_violation(scan, "|a| above a_max", t, a);
    }
}

template <typename Scalar>
void record_jerk(envelope_scan<Scalar>& scan, Scalar t_prev, Scalar t, Scalar a_prev, Scalar a)
{
    auto const eps = std::numeric_limits<Scalar>::epsilon();
    auto const& limits = scan.limits;
    auto const time_scale = std::max(std::abs(t), scan.planner.diagnostics().planned_duration);
    auto const allowance = Scalar{2}
                           * (static_cast<Scalar>(time_rounding_ops) * eps * limits.j_max * time_scale
                              + static_cast<Scalar>(acceleration_rounding_ops) * eps * limits.a_max);
    scan.largest_allowance = std::max(scan.largest_allowance, allowance);
    auto const bound = limits.j_max * (t - t_prev) + allowance;
    if (std::abs(a - a_prev) > bound) {
        ++scan.counts.jerk_violations;
        report_violation(scan, "|a_k - a_(k-1)| above j_max dt plus the allowance", t, std::abs(a - a_prev) - bound);
    }
}

template <typename Scalar>
void scan_dense(envelope_scan<Scalar>& scan, Scalar t_begin)
{
    auto const dt = static_cast<Scalar>(grid_step);
    auto const steps = static_cast<std::int32_t>(std::lround(dense_seconds / grid_step));
    Scalar t_prev{};
    Scalar a_prev{};
    for (std::int32_t k = 0; k <= steps; ++k) {
        auto const t = t_begin + static_cast<Scalar>(k) * dt;
        auto const point = scan.planner.sample(t);
        record_sample(scan, t, point);
        if (k > 0) {
            record_jerk(scan, t_prev, t, a_prev, point.acceleration[0]);
        }
        t_prev = t;
        a_prev = point.acceleration[0];
    }
}

// Dense over the first and last stretch of the plan, where every phase but the cruise lives, and every second between.
template <typename Scalar>
void scan_plan(envelope_scan<Scalar>& scan, Scalar t0)
{
    auto const duration = scan.planner.diagnostics().planned_duration;
    auto const dense = static_cast<Scalar>(dense_seconds);
    scan_dense(scan, t0);
    for (std::int32_t k = dense_seconds; static_cast<Scalar>(k) + dense < duration; ++k) {
        auto const t = t0 + static_cast<Scalar>(k);
        record_sample(scan, t, scan.planner.sample(t));
    }
    scan_dense(scan, t0 + std::max(duration - dense, Scalar{0}));
}

}

#endif
