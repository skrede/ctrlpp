#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_RETARGET_WINDOWS_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_RETARGET_WINDOWS_H

#include "common.h"

#include <cmath>
#include <cstdint>

namespace planner_overshoot_sweep
{

// The last three are sets on which a cruise integrated forward through its ramp rounds above v_max.
inline constexpr limit_set retarget_limits[] = {
    {1.0, 1.0, 1.0},       {5.0, 10.0, 50.0}, {2.0, 3.0, 7.0},     {0.7, 1.3, 2.9},
    {10.0, 100.0, 1000.0}, {1.0, 5.0, 50.0},  {40.0, 20.0, 100.0}, {1e-3, 1e-2, 1e-1},
};

inline constexpr double grid_step = 1e-3;
inline constexpr std::int32_t instants_per_window = 40;
inline constexpr std::int32_t dense_seconds = 20;

// Elapsed time reaches a sample through at most 13 roundings, and the sampled acceleration carries at most 4
// since the last pinned phase end; a finite-difference jerk is admitted twice the excess those counts allow.
inline constexpr std::int32_t time_rounding_ops = 13;
inline constexpr std::int32_t acceleration_rounding_ops = 4;

enum class retarget_window : std::uint8_t
{
    cruise,
    acceleration_ramp,
    acceleration_nulling,
};

inline constexpr retarget_window retarget_windows[] = {retarget_window::cruise, retarget_window::acceleration_ramp,
                                                       retarget_window::acceleration_nulling};

inline auto window_name(retarget_window window) -> char const*
{
    switch (window) {
    case retarget_window::cruise:               return "cruise";
    case retarget_window::acceleration_ramp:    return "acceleration ramp";
    case retarget_window::acceleration_nulling: return "acceleration nulling";
    }
    return "unnamed";
}

template <typename Scalar>
struct window_span
{
    Scalar begin;
    Scalar end;
};

// @cite biagiotti2009 -- Sec. 3.4.3, rest-to-rest ramp with and without a_max reached
template <typename Scalar>
auto window_span_of(retarget_window window, typed_limits<Scalar> const& limits, Scalar dt) -> window_span<Scalar>
{
    auto T_j = limits.a_max / limits.j_max;
    auto T_a = T_j + limits.v_max / limits.a_max;
    if (limits.v_max * limits.j_max < limits.a_max * limits.a_max) {
        T_j = std::sqrt(limits.v_max / limits.j_max);
        T_a = Scalar{2} * T_j;
    }
    switch (window) {
    case retarget_window::cruise:               return {Scalar{1.5} * T_a, Scalar{1.5} * T_a + Scalar{40} * dt};
    case retarget_window::acceleration_ramp:    return {dt, T_a - dt};
    case retarget_window::acceleration_nulling: break;
    }
    return {T_a - T_j, T_a - dt};
}

template <typename Scalar>
auto window_instant(window_span<Scalar> const& span, Scalar dt, std::int32_t i, std::int32_t count) -> Scalar
{
    auto const first = static_cast<std::int32_t>(std::ceil(span.begin / dt));
    auto const last = static_cast<std::int32_t>(std::floor(span.end / dt));
    auto const k = first + (last - first) * i / (count - 1);
    return static_cast<Scalar>(k) * dt;
}

}

#endif
