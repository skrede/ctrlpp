#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_SCAN_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_SCAN_H

#include "common.h"
#include "envelope_legs.h"

#include "ctrlpp/trajectory/trajectory_types.h"

#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdint>
#include <algorithm>

namespace planner_overshoot_sweep
{

inline constexpr std::int32_t switch_window_steps = 64;

struct switch_counts
{
    std::int64_t configurations;
    std::int64_t not_measured;
    std::int64_t instants;
    std::int64_t windows_missing;
    std::int64_t windows_not_located;
    std::int64_t velocity_violations;
    std::int64_t acceleration_violations;
};

enum class switch_evidence : std::uint8_t
{
    none,
    acceleration_changes,
    limit_then_zero,
};

template <typename Scalar>
struct switch_scan
{
    switch_counts& counts;
    char const* configuration;
    Scalar v_max;
    Scalar a_bound;
    switch_evidence evidence;
    bool closed_form;
    bool print_points;
};

template <typename Scalar>
struct instant_span
{
    Scalar lo;
    Scalar hi;
};

struct window_trace
{
    bool acceleration_changed;
    bool limit_reached;
    bool limit_then_zero;
};

template <typename Scalar>
void report_switch_excess(switch_scan<Scalar> const& scan, char const* what, Scalar t, Scalar value, Scalar limit)
{
    if (scan.print_points) {
        std::printf("      SWITCH  [%s]  %s  t=%a  value=%.17g  limit=%.17g  excess %lld ulp\n", scan.configuration, what,
                    static_cast<double>(t), static_cast<double>(value), static_cast<double>(limit),
                    static_cast<long long>(ulps_between(limit, value)));
    }
}

template <typename Scalar>
void record_switch_sample(switch_scan<Scalar>& scan, Scalar t, trajectory_point<Scalar, 1> const& point)
{
    ++scan.counts.instants;
    auto const v = std::abs(point.velocity[0]);
    auto const a = std::abs(point.acceleration[0]);
    if (v > scan.v_max) {
        ++scan.counts.velocity_violations;
        report_switch_excess(scan, "|v|", t, v, scan.v_max);
    }
    if (a > scan.a_bound) {
        ++scan.counts.acceleration_violations;
        report_switch_excess(scan, "|a|", t, a, scan.a_bound);
    }
}

template <typename Scalar>
void trace_acceleration(window_trace& trace, Scalar first, Scalar a, Scalar a_bound)
{
    trace.acceleration_changed = trace.acceleration_changed || a != first;
    trace.limit_then_zero = trace.limit_then_zero || (trace.limit_reached && a == Scalar{0});
    trace.limit_reached = trace.limit_reached || std::abs(a) == a_bound;
}

inline auto witnessed(switch_evidence evidence, window_trace const& trace) -> bool
{
    switch (evidence) {
    case switch_evidence::none:                 return true;
    case switch_evidence::acceleration_changes: return trace.acceleration_changed;
    case switch_evidence::limit_then_zero:      return trace.limit_then_zero;
    }
    return false;
}

template <typename Scalar, typename Sample>
void scan_span(switch_scan<Scalar>& scan, instant_span<Scalar> const& span, Sample const& sample)
{
    window_trace trace{};
    auto const first = sample(span.lo).acceleration[0];
    for (auto t = span.lo;; t = std::nextafter(t, span.hi)) {
        auto const point = sample(t);
        record_switch_sample(scan, t, point);
        trace_acceleration(trace, first, point.acceleration[0], scan.a_bound);
        if (!(t < span.hi)) {
            break;
        }
    }
    if (witnessed(scan.evidence, trace)) {
        return;
    }
    // A window placed by a closed form the planner does not mirror exactly is reported, never a zero claim.
    ++(scan.closed_form ? scan.counts.windows_not_located : scan.counts.windows_missing);
    if (scan.print_points) {
        std::printf("      %s  [%s]  window %a .. %a\n", scan.closed_form ? "NOT LOCATED" : "MISSING", scan.configuration,
                    static_cast<double>(span.lo), static_cast<double>(span.hi));
    }
}

template <typename Scalar>
auto span_around(Scalar center, Scalar lowest, Scalar highest) -> instant_span<Scalar>
{
    auto lo = std::clamp(center, lowest, highest);
    auto hi = lo;
    for (std::int32_t i = 0; i < switch_window_steps; ++i) {
        lo = std::nextafter(lo, lowest);
        hi = std::nextafter(hi, highest);
    }
    return {lo, hi};
}

// Windows that overlap are merged, so every instant is sampled once and each merged window holds at least one switch.
template <typename Scalar, typename Sample>
void scan_switches(switch_scan<Scalar>& scan, std::vector<Scalar> centers, instant_span<Scalar> const& domain, Sample const& sample)
{
    std::sort(centers.begin(), centers.end());
    auto open = span_around(centers.front(), domain.lo, domain.hi);
    for (auto const center : centers) {
        auto const span = span_around(center, domain.lo, domain.hi);
        if (span.lo <= open.hi) {
            open.hi = std::max(open.hi, span.hi);
            continue;
        }
        scan_span(scan, open, sample);
        open = span;
    }
    scan_span(scan, open, sample);
}

inline void report_switch_block(char const* name, char const* type, std::uint32_t seed, switch_counts const& counts)
{
    std::printf("\n  Phase switches, %s, %s totals\n", name, type);
    report_count("switch seed", seed);
    report_count("switch configurations", counts.configurations);
    report_count("switch configurations not measured", counts.not_measured);
    report_count("switch instants sampled", counts.instants);
    report_count("switch windows missing the switch", counts.windows_missing);
    report_count("switch windows not located", counts.windows_not_located);
    report_count("switch velocity samples above the limit", counts.velocity_violations);
    report_count("switch acceleration samples above the limit", counts.acceleration_violations);
}

}

#endif
