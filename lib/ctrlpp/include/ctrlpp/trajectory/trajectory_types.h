#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_TYPES_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_TYPES_H

/// @brief Core output and error types for trajectory generation.
///
/// trajectory_point<Scalar, ND> holds position, velocity, and acceleration vectors
/// for an ND-dimensional trajectory. path_point<Scalar> holds normalized scalar
/// values (q, dq, ddq, dddq) for a path mapping tau in [0,1]. spline_error and
/// trajectory_error enumerate the structured failure modes of the spline and
/// trajectory factories; each forms the error channel of a
/// `ctrlpp::expected<T, E>` contract.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 5.2.1, eq. (2.16)

#include "ctrlpp/types.h"

#include <cstddef>

namespace ctrlpp
{

/// @brief Structured failure modes for the spline factories: `cubic_spline`,
/// `smoothing_spline`, `bspline_trajectory`, and `make_bspline_interpolation`.
///
///  * too_few_points             : fewer waypoints than the factory minimum
///                                 (2 for cubic and smoothing splines,
///                                 Degree + 1 for B-spline interpolation).
///  * size_mismatch              : times and positions differ in length.
///  * non_increasing_times       : knot times are not strictly increasing.
///  * periodic_endpoint_mismatch : periodic boundary conditions require the first
///                                 and last positions to match within the endpoint
///                                 rounding budget.
///  * periodic_too_few_points    : periodic boundary conditions require at least
///                                 3 waypoints.
///  * too_few_control_points     : a B-spline of degree p requires at least p + 1
///                                 control points.
///  * bad_knot_count             : the knot vector size differs from
///                                 control_points.size() + Degree + 1.
///  * non_monotonic_knots        : the knot vector is not non-decreasing.
///  * mu_out_of_range            : smoothing parameter mu lies outside (0, 1].
enum class spline_error
{
    too_few_points,
    size_mismatch,
    non_increasing_times,
    periodic_endpoint_mismatch,
    periodic_too_few_points,
    too_few_control_points,
    bad_knot_count,
    non_monotonic_knots,
    mu_out_of_range,
};

/// @brief Structured failure modes for the point-to-point trajectory factories
/// and the online trajectory planners.
///
///  * non_positive_velocity_limit     : the velocity limit must be positive.
///  * non_positive_acceleration_limit : the acceleration limit must be positive.
///  * non_positive_jerk_limit         : the jerk limit must be positive.
///  * non_positive_duration           : the requested duration must be positive.
///  * non_finite_input                : a boundary value or limit is NaN/Inf.
///  * boundary_velocity_exceeds_limit : a commanded boundary velocity is larger
///                                      in magnitude than the velocity limit it
///                                      is commanded under. The limit is a
///                                      precondition of the point-to-point
///                                      profiles, not a value they raise to fit:
///                                      raising it would violate a bound the
///                                      caller asked for, and honoring it would
///                                      require a ramp that runs backwards in
///                                      time.
///  * unrepresentable_duration        : a duration is not representable in the
///                                      scalar type. Three cases, all of them
///                                      facts about the arithmetic rather than
///                                      about the kinematics: a duration of the
///                                      constructed profile left the finite
///                                      range; the total underflowed to zero on
///                                      a command with a nonzero displacement,
///                                      which would report an instantaneous
///                                      traversal; or a requested retiming
///                                      differs from a duration the profile
///                                      already realizes by less than the
///                                      rounding of the expression that would
///                                      solve for it, so the request cannot be
///                                      told apart from that duration. The third
///                                      case says a profile may well exist and
///                                      the arithmetic cannot locate it, which
///                                      is why it is reported separately from
///                                      unreachable_duration: the caller should
///                                      change the request, not the limits.
///  * unreachable_boundary_velocity   : the commanded displacement is smaller
///                                      than the distance the fastest admissible
///                                      transition between the two boundary
///                                      velocities already sweeps, so no profile
///                                      of the requested shape realizes it.
///  * duration_shorter_than_current   : time rescaling only slows a profile down.
///                                      The profile already runs at the fastest
///                                      shape its limits allow, so a duration
///                                      below the current one is not realizable.
///  * unreachable_duration            : the requested duration lies outside the
///                                      set the commanded displacement, the
///                                      kinematic limits, and the boundary
///                                      velocities can realize together.
enum class trajectory_error
{
    non_positive_velocity_limit,
    non_positive_acceleration_limit,
    non_positive_jerk_limit,
    non_positive_duration,
    non_finite_input,
    boundary_velocity_exceeds_limit,
    unrepresentable_duration,
    unreachable_boundary_velocity,
    duration_shorter_than_current,
    unreachable_duration,
};

/// @brief Point on an ND-dimensional trajectory with position, velocity, acceleration.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Sec. 2.1, eq. (2.1)-(2.3) -- trajectory representation
///   with position q(t), velocity q'(t), acceleration q''(t)
template <typename Scalar, std::size_t ND>
struct trajectory_point
{
    Vector<Scalar, ND> position{};
    Vector<Scalar, ND> velocity{};
    Vector<Scalar, ND> acceleration{};
};

/// @brief Normalized path output for tau in [0,1].
///
/// Fields represent normalized position, velocity, acceleration, and jerk.
/// Physical values are obtained via scaling: vel = h/T * dq, acc = h/T^2 * ddq, etc.
template <typename Scalar>
struct path_point
{
    Scalar q{};
    Scalar dq{};
    Scalar ddq{};
    Scalar dddq{};
};

}

#endif
