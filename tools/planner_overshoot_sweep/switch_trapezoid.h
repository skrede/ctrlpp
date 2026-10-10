#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_TRAPEZOID_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_TRAPEZOID_H

#include "common.h"
#include "switch_scan.h"
#include "switch_shapes.h"
#include "switch_population.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <cstdio>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <type_traits>

namespace planner_overshoot_sweep
{

inline constexpr double rescale_factors[] = {1.5, 3.0};

// Fractions of the ground the transition between the two boundary velocities sweeps. Exact powers of two, so
// each command is realizable only at an acceleration above the commanded one, by about the reciprocal.
inline constexpr double raising_fractions[] = {0.5, 0x1p-20};

struct policy_counts
{
    std::int64_t refused;
    std::int64_t disagreements;
};

template <typename Scalar>
using trapezoid_result = ctrlpp::expected<trapezoidal_trajectory<Scalar>, trajectory_error>;

template <typename Scalar>
auto create_trapezoid(shape_command<Scalar> const& c, ctrlpp::acceleration_limit_policy policy) -> trapezoid_result<Scalar>
{
    return trapezoidal_trajectory<Scalar>::create({.q0 = c.q0, .q1 = c.q1, .v_max = c.limits.v_max, .a_max = c.limits.a_max,
                                                   .v0 = c.v0, .v1 = c.v1, .acceleration_policy = policy});
}

template <typename Scalar>
auto same_profile(trapezoidal_trajectory<Scalar> const& lhs, trapezoidal_trajectory<Scalar> const& rhs) -> bool
{
    return lhs.phase_durations() == rhs.phase_durations() && lhs.duration() == rhs.duration()
           && lhs.peak_velocity() == rhs.peak_velocity()
           && lhs.disposition().realized_acceleration == rhs.disposition().realized_acceleration;
}

// The default refuses exactly where raise_to_feasible builds above the commanded acceleration; every other
// outcome of the two policies is identical.
template <typename Scalar>
auto policies_agree(trapezoid_result<Scalar> const& raised, trapezoid_result<Scalar> const& refused) -> bool
{
    bool const above = raised.has_value()
                       && raised->disposition().realized_acceleration != raised->disposition().commanded_acceleration;
    if (above) {
        return !refused.has_value() && refused.error() == trajectory_error::acceleration_limit_infeasible;
    }
    if (raised.has_value() != refused.has_value()) {
        return false;
    }
    return raised.has_value() ? same_profile(*raised, *refused) : raised.error() == refused.error();
}

template <typename Scalar>
void measure_trapezoid(shape_command<Scalar> const& c, char const* tag, shape_leg& leg, policy_counts& policies)
{
    ++leg.counts.configurations;
    char configuration[320];
    describe(c, tag, configuration);
    auto const created = create_trapezoid(c, ctrlpp::acceleration_limit_policy::raise_to_feasible);
    auto const refused = create_trapezoid(c, ctrlpp::acceleration_limit_policy::refuse);
    policies.refused += (!refused.has_value() && refused.error() == trajectory_error::acceleration_limit_infeasible) ? 1 : 0;
    policies.disagreements += policies_agree(created, refused) ? 0 : 1;
    if (!created.has_value()) {
        ++leg.counts.not_measured;
        return;
    }
    scan_shape(*created, c.limits, configuration, leg);
    for (auto const factor : rescale_factors) {
        auto rescaled = *created;
        auto const T_new = static_cast<Scalar>(factor) * created->duration();
        if (rescaled.can_rescale_to(T_new).has_value() && rescaled.rescale_to(T_new).has_value()) {
            ++leg.rescaled;
            char stretched[352];
            std::snprintf(stretched, sizeof(stretched), "%s rescaled x%g", configuration, factor);
            scan_shape(rescaled, c.limits, stretched, leg);
        }
    }
}

// In the positive frame the two ramps sweep at least (v_hi - v_lo)(v_hi + v_lo) / (2 a), B&M eq. (3.14).
template <typename Scalar>
auto raising_variants(shape_command<Scalar> const& c) -> std::vector<shape_command<Scalar>>
{
    auto const sigma = c.q1 < c.q0 ? Scalar{-1} : Scalar{1};
    auto const v_hi = std::max(sigma * c.v0, sigma * c.v1);
    auto const v_lo = std::min(sigma * c.v0, sigma * c.v1);
    auto const ground = (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * c.limits.a_max);
    std::vector<shape_command<Scalar>> variants;
    if (ground > Scalar{0}) {
        for (auto const fraction : raising_fractions) {
            variants.push_back({c.q0, c.q0 + sigma * static_cast<Scalar>(fraction) * ground, c.v0, c.v1, c.limits});
        }
    }
    return variants;
}

// Configurations on which an earlier random search sampled a value above a limit next to a phase switch.
template <typename Scalar>
auto trapezoid_reproducers() -> std::vector<shape_command<Scalar>>
{
    if constexpr (std::is_same_v<Scalar, double>) {
        return {{0.35790432359834429, 9621.0167335908081, 0.0028366887302119395, 0.006388618899210774,
                 {0.006846007648171031, 0.20062093968092909, 0.0}}};
    } else {
        return {};
    }
}

template <typename Scalar>
void run_trapezoid_switch_leg(bool print_points)
{
    auto const seed = switch_seed<Scalar>(switch_leg::trapezoidal);
    fraction_source<Scalar> source{seed};
    shape_leg leg{{}, 0, print_points};
    policy_counts policies{};
    std::int64_t raising = 0;
    for_each_trapezoid_command(source, [&](shape_command<Scalar> const& c) {
        measure_trapezoid(c, "", leg, policies);
        for (auto const& variant : raising_variants(c)) {
            ++raising;
            measure_trapezoid(variant, "raising ", leg, policies);
        }
    });
    for (auto const& c : trapezoid_reproducers<Scalar>()) {
        measure_trapezoid(c, "reproducer ", leg, policies);
    }
    if (print_points) {
        char parts[200];
        std::snprintf(parts, sizeof(parts), "commands %lld of which raising %lld, rescaled profiles %lld, reproducers %zu",
                      static_cast<long long>(leg.counts.configurations), static_cast<long long>(raising),
                      static_cast<long long>(leg.rescaled), trapezoid_reproducers<Scalar>().size());
        report_population("trapezoidal", scalar_name<Scalar>(), parts);
    }
    report_switch_block("trapezoidal", scalar_name<Scalar>(), seed, leg.counts);
    report_count("refusals under the default policy", policies.refused);
    report_count("default outcomes unlike raise_to_feasible", policies.disagreements);
}

}

#endif
