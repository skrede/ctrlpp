#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_PRECONDITION_BUILD_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_PRECONDITION_BUILD_H

#include "common.h"
#include "retarget_windows.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <type_traits>

namespace planner_overshoot_sweep
{

inline constexpr std::int32_t lowest_speed_step = -3;
inline constexpr std::int32_t highest_speed_step = 6;
inline constexpr std::uint32_t knob_values[] = {0, 1, 2, 4};
inline constexpr auto error_count = static_cast<std::size_t>(trajectory_error::unreachable_duration) + 1;

struct precondition_counts
{
    std::int64_t commands;
    std::int64_t builds;
    std::array<std::int64_t, error_count> refusals;
    std::int64_t peak_above_limit;
    std::int64_t refused_at_or_below;
    std::int64_t refused_within_allowance;
    std::int64_t admitted_beyond_allowance;
    std::int64_t clamped_exact;
    std::int64_t clamped_off;
    double largest_peak_ratio;
};

inline void accumulate(precondition_counts& total, precondition_counts const& row)
{
    total.commands += row.commands;
    total.builds += row.builds;
    for (std::size_t e = 0; e < error_count; ++e) {
        total.refusals[e] += row.refusals[e];
    }
    total.peak_above_limit += row.peak_above_limit;
    total.refused_at_or_below += row.refused_at_or_below;
    total.refused_within_allowance += row.refused_within_allowance;
    total.admitted_beyond_allowance += row.admitted_beyond_allowance;
    total.clamped_exact += row.clamped_exact;
    total.clamped_off += row.clamped_off;
    total.largest_peak_ratio = std::max(total.largest_peak_ratio, row.largest_peak_ratio);
}

template <typename Scalar>
struct boundary_command
{
    std::int32_t step;
    Scalar speed;
    Scalar sigma;
    bool final_role;
    std::uint32_t knob;
};

template <typename Scalar>
auto stepped_speed(Scalar v_max, std::int32_t step) -> Scalar
{
    auto const toward = step < 0 ? Scalar{0} : std::numeric_limits<Scalar>::infinity();
    auto speed = v_max;
    for (std::int32_t k = 0; k < std::abs(step); ++k) {
        speed = std::nextafter(speed, toward);
    }
    return speed;
}

template <typename Scalar>
auto allowance_of(boundary_command<Scalar> const& command, Scalar v_max) -> Scalar
{
    return static_cast<Scalar>(command.knob) * std::numeric_limits<Scalar>::epsilon() * v_max;
}

template <typename Shape, typename Scalar>
auto build_shape(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command)
    -> ctrlpp::expected<Shape, trajectory_error>
{
    auto const v = command.sigma * command.speed;
    auto const q1 = command.sigma * Scalar{10} * limits.v_max;
    auto const v0 = command.final_role ? Scalar{0} : v;
    auto const v1 = command.final_role ? v : Scalar{0};
    if constexpr (std::is_same_v<Shape, double_s_trajectory<Scalar>>) {
        return Shape::create({.q0 = Scalar{0}, .q1 = q1, .v_max = limits.v_max, .a_max = limits.a_max,
                              .j_max = limits.j_max, .v0 = v0, .v1 = v1,
                              .boundary_velocity_rounding_ops = command.knob});
    } else {
        return Shape::create({.q0 = Scalar{0}, .q1 = q1, .v_max = limits.v_max, .a_max = limits.a_max,
                              .v0 = v0, .v1 = v1, .boundary_velocity_rounding_ops = command.knob});
    }
}

template <typename Scalar>
auto realized_peak(double_s_trajectory<Scalar> const& shape) -> Scalar
{
    return std::abs(shape.peak_velocity());
}

template <typename Scalar>
auto realized_peak(trapezoidal_trajectory<Scalar> const& shape) -> Scalar
{
    auto const dt = static_cast<Scalar>(grid_step);
    auto peak = std::abs(shape.evaluate(shape.duration()).velocity[0]);
    for (std::int32_t k = 0; static_cast<Scalar>(k) * dt < shape.duration(); ++k) {
        peak = std::max(peak, std::abs(shape.evaluate(static_cast<Scalar>(k) * dt).velocity[0]));
    }
    return peak;
}

template <typename Scalar>
void print_command(char const* configuration, boundary_command<Scalar> const& command, char const* outcome)
{
    std::printf("      [%s  speed %+d ulp  sigma %+d  %s  knob %u]  %s\n", configuration, command.step,
                static_cast<std::int32_t>(command.sigma), command.final_role ? "v1" : "v0", command.knob, outcome);
}

template <typename Scalar>
void record_refusal(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command,
                    trajectory_error error, precondition_counts& counts)
{
    ++counts.refusals[static_cast<std::size_t>(error)];
    if (command.speed <= limits.v_max) {
        ++counts.refused_at_or_below;
    } else if (command.speed - limits.v_max <= allowance_of(command, limits.v_max)) {
        ++counts.refused_within_allowance;
    }
}

template <typename Shape, typename Scalar>
auto record_build(typed_limits<Scalar> const& limits, boundary_command<Scalar> const& command, Shape const& shape,
                  precondition_counts& counts) -> Scalar
{
    ++counts.builds;
    auto const peak = realized_peak(shape);
    counts.peak_above_limit += (peak > limits.v_max) ? 1 : 0;
    counts.largest_peak_ratio = std::max(counts.largest_peak_ratio, static_cast<double>(peak / limits.v_max));
    if (command.speed > limits.v_max) {
        counts.admitted_beyond_allowance += (command.speed - limits.v_max > allowance_of(command, limits.v_max)) ? 1 : 0;
        auto const end = shape.evaluate(command.final_role ? shape.duration() : Scalar{0}).velocity[0];
        ++(end == command.sigma * limits.v_max ? counts.clamped_exact : counts.clamped_off);
    }
    return peak;
}

}

#endif
