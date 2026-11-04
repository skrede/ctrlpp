#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_ERROR_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_ERROR_H

namespace ctrlpp
{

enum class spline_error
{
    too_few_points,
    size_mismatch,
    non_finite_input,
    non_increasing_times,
    periodic_endpoint_mismatch,
    periodic_too_few_points,
    too_few_control_points,
    bad_knot_count,
    non_monotonic_knots,
    mu_out_of_range,
    unrepresentable_spline,
};

enum class trajectory_error
{
    non_positive_velocity_limit,
    non_positive_acceleration_limit,
    non_positive_jerk_limit,
    non_positive_duration,
    non_finite_input,
    boundary_velocity_exceeds_limit,
    boundary_velocity_rounding_ops_out_of_range,
    unrepresentable_duration,
    unreachable_boundary_velocity,
    acceleration_limit_infeasible,
    duration_shorter_than_current,
    unreachable_duration,
};

}

#endif
