#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_SHAPES_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_SHAPES_H

#include "common.h"
#include "switch_scan.h"
#include "switch_population.h"

#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <cstdio>
#include <vector>
#include <cstdint>
#include <type_traits>

namespace planner_overshoot_sweep
{

struct shape_leg
{
    switch_counts counts;
    std::int64_t rescaled;
    bool print_points;
};

template <typename Scalar>
void describe(shape_command<Scalar> const& c, char const* tag, char (&out)[320])
{
    std::snprintf(out, sizeof(out), "%sq0=%a q1=%a v=%a a=%a j=%a v0=%a v1=%a", tag, static_cast<double>(c.q0),
                  static_cast<double>(c.q1), static_cast<double>(c.limits.v_max), static_cast<double>(c.limits.a_max),
                  static_cast<double>(c.limits.j_max), static_cast<double>(c.v0), static_cast<double>(c.v1));
}

template <typename Scalar>
auto switch_centers(trapezoidal_trajectory<Scalar> const& shape) -> std::vector<Scalar>
{
    auto const p = shape.phase_durations();
    return {p[0], p[0] + p[1], p[0] + p[1] + p[2], shape.duration()};
}

template <typename Scalar>
auto switch_centers(double_s_trajectory<Scalar> const& shape) -> std::vector<Scalar>
{
    std::vector<Scalar> centers{shape.duration()};
    auto cumulative = Scalar{0};
    for (auto const phase : shape.phase_durations()) {
        cumulative += phase;
        centers.push_back(cumulative);
    }
    return centers;
}

template <typename Scalar>
auto acceleration_bound(trapezoidal_trajectory<Scalar> const& shape, typed_limits<Scalar> const&) -> Scalar
{
    return shape.disposition().realized_acceleration;
}

template <typename Scalar>
auto acceleration_bound(double_s_trajectory<Scalar> const&, typed_limits<Scalar> const& limits) -> Scalar
{
    return limits.a_max;
}

template <typename Shape, typename Scalar>
void scan_shape(Shape const& shape, typed_limits<Scalar> const& limits, char const* configuration, shape_leg& leg)
{
    switch_scan<Scalar> scan{leg.counts, configuration,        limits.v_max, acceleration_bound(shape, limits),
                             switch_evidence::none, false,          leg.print_points};
    auto const sample = [&](Scalar t) { return shape.evaluate(t); };
    scan_switches(scan, switch_centers(shape), instant_span<Scalar>{Scalar{0}, shape.duration()}, sample);
}

template <typename Scalar>
void measure_double_s(shape_command<Scalar> const& c, char const* tag, shape_leg& leg)
{
    ++leg.counts.configurations;
    char configuration[320];
    describe(c, tag, configuration);
    auto const created = double_s_trajectory<Scalar>::create({.q0 = c.q0, .q1 = c.q1, .v_max = c.limits.v_max, .a_max = c.limits.a_max,
                                                              .j_max = c.limits.j_max, .v0 = c.v0, .v1 = c.v1});
    if (!created.has_value()) {
        ++leg.counts.not_measured;
        return;
    }
    scan_shape(*created, c.limits, configuration, leg);
}

template <typename Scalar>
auto double_s_reproducers() -> std::vector<shape_command<Scalar>>
{
    if constexpr (std::is_same_v<Scalar, double>) {
        auto const v = 0.086115804201940746;
        auto const d = 0.50507005729686805;
        return {{0.0, 1015.7326884443336, 0.0, -v, {v, 0.0037083513761102535, 4124533871.0118299}},
                {0.0, 1000.0, d, 0.0, {2.0 * d, 0.029259168123834691, 0.0016950102405212067}}};
    } else {
        auto const v = static_cast<float>(0.0044597648084163666);
        return {{static_cast<float>(0.00033788720611482859), static_cast<float>(-5.8095717430114746),
                 static_cast<float>(-0.0035865840036422014), v,
                 {v, static_cast<float>(0.0015443010488525033), static_cast<float>(30.257682800292969)}}};
    }
}

inline void report_population(char const* name, char const* type, char const* parts)
{
    std::printf("  [%s, %s population]  %s\n", name, type, parts);
}

template <typename Scalar>
void run_double_s_switch_leg(bool print_points)
{
    auto const seed = switch_seed<Scalar>(switch_leg::double_s);
    fraction_source<Scalar> source{seed};
    shape_leg leg{{}, 0, print_points};
    auto const measure = [&](shape_command<Scalar> const& c) { measure_double_s(c, "", leg); };
    for_each_seeded_double_s_command(source, measure);
    auto const seeded = leg.counts.configurations;
    for_each_long_double_s_command(source, measure);
    auto const long_profiles = leg.counts.configurations - seeded;
    for_each_threshold_double_s_command<Scalar>(measure);
    auto const threshold = leg.counts.configurations - seeded - long_profiles;
    for (auto const& c : double_s_reproducers<Scalar>()) {
        measure_double_s(c, "reproducer ", leg);
    }
    if (print_points) {
        auto const m = long_exponents<Scalar>();
        char parts[200];
        std::snprintf(parts, sizeof(parts), "seeded %lld, long %lld at 2^m for m = %d %d %d, threshold %lld, reproducers %zu",
                      static_cast<long long>(seeded), static_cast<long long>(long_profiles), m[0], m[1], m[2],
                      static_cast<long long>(threshold), double_s_reproducers<Scalar>().size());
        report_population("double_s", scalar_name<Scalar>(), parts);
    }
    report_switch_block("double_s", scalar_name<Scalar>(), seed, leg.counts);
}

}

#endif
