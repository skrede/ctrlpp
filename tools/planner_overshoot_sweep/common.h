#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_COMMON_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_COMMON_H

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <bit>
#include <cstdint>
#include <type_traits>

namespace planner_overshoot_sweep
{

using ctrlpp::double_s_trajectory;
using ctrlpp::online_planner_2nd;
using ctrlpp::online_planner_3rd;
using ctrlpp::online_planner_substitution_reason;
using ctrlpp::trajectory_error;
using ctrlpp::trajectory_point;
using ctrlpp::trapezoidal_trajectory;

/// How far past each located boundary the exhaustive walk extends, in
/// representable targets. Wide enough that a boundary landing a few values away
/// from where the bisection put it is still inside the walk, and narrow enough
/// that the walk's own count stays readable.
inline constexpr std::int32_t walk_margin_ulps = 64;

struct limit_set
{
    double v_max{};
    double a_max{};
    double j_max{};
};

template <typename Scalar>
struct typed_limits
{
    Scalar v_max;
    Scalar a_max;
    Scalar j_max;
};

template <typename Scalar>
auto typed(limit_set const& limits) -> typed_limits<Scalar>
{
    return {static_cast<Scalar>(limits.v_max), static_cast<Scalar>(limits.a_max),
            static_cast<Scalar>(limits.j_max)};
}

template <typename Scalar>
auto scalar_name() -> char const*
{
    return std::is_same_v<Scalar, float> ? "float" : "double";
}

template <typename Scalar>
using target_index_t = std::conditional_t<sizeof(Scalar) == 4, std::int32_t, std::int64_t>;

template <typename Scalar>
using target_bits_t = std::conditional_t<sizeof(Scalar) == 4, std::uint32_t, std::uint64_t>;

/// Position of a strictly positive value in the ordered sequence of the values
/// its type represents. Adjacent representable values differ by one here, so a
/// difference of two such positions is a count of units in the last place
/// rather than a length that has to be divided by one.
template <typename Scalar>
auto target_index(Scalar value) -> target_index_t<Scalar>
{
    return static_cast<target_index_t<Scalar>>(std::bit_cast<target_bits_t<Scalar>>(value));
}

template <typename Scalar>
auto target_at(target_index_t<Scalar> index) -> Scalar
{
    return std::bit_cast<Scalar>(static_cast<target_bits_t<Scalar>>(index));
}

template <typename Scalar>
auto ulps_between(Scalar low, Scalar high) -> target_index_t<Scalar>
{
    return target_index(high) - target_index(low);
}

inline auto error_name(trajectory_error error) -> char const*
{
    switch (error) {
    case trajectory_error::non_positive_velocity_limit:     return "non_positive_velocity_limit";
    case trajectory_error::non_positive_acceleration_limit: return "non_positive_acceleration_limit";
    case trajectory_error::non_positive_jerk_limit:         return "non_positive_jerk_limit";
    case trajectory_error::non_positive_duration:           return "non_positive_duration";
    case trajectory_error::non_finite_input:                return "non_finite_input";
    case trajectory_error::boundary_velocity_exceeds_limit: return "boundary_velocity_exceeds_limit";
    case trajectory_error::boundary_velocity_rounding_ops_out_of_range:
        return "boundary_velocity_rounding_ops_out_of_range";
    case trajectory_error::unrepresentable_duration:        return "unrepresentable_duration";
    case trajectory_error::unreachable_boundary_velocity:   return "unreachable_boundary_velocity";
    case trajectory_error::duration_shorter_than_current:   return "duration_shorter_than_current";
    case trajectory_error::unreachable_duration:            return "unreachable_duration";
    }
    return "unnamed";
}

inline auto reason_name(online_planner_substitution_reason reason) -> char const*
{
    switch (reason) {
    case online_planner_substitution_reason::none:                  return "none";
    case online_planner_substitution_reason::reversal_or_overshoot: return "reversal_or_overshoot";
    case online_planner_substitution_reason::carry_velocity_shape_unavailable:
        return "carry_velocity_shape_unavailable";
    }
    return "unnamed";
}

}

#endif
