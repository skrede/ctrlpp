#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_ERROR_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_ERROR_H

/// @brief Error types for trajectory generation.
///
/// spline_error and trajectory_error enumerate the structured failure modes of the
/// spline and trajectory factories; each forms the error channel of a
/// `ctrlpp::expected<T, E>` contract.

namespace ctrlpp
{

/// @brief Structured failure modes for the spline factories: `cubic_spline`,
/// `smoothing_spline`, `bspline_trajectory`, and `make_bspline_interpolation`.
///
///  * too_few_points             : fewer waypoints than the factory minimum
///                                 (2 for cubic and smoothing splines,
///                                 Degree + 1 for B-spline interpolation).
///  * size_mismatch              : times and positions differ in length.
///  * non_finite_input           : a time, position, boundary velocity, control
///                                 point, or knot is NaN or infinite.
///  * non_increasing_times       : knot times are not strictly increasing, or a
///                                 B-spline's active domain has zero length.
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
///  * unrepresentable_spline     : the configuration is well posed, and the
///                                 arithmetic that would realize it is not
///                                 representable in Scalar. It is reported apart
///                                 from every rejection above because it is a
///                                 statement about the type rather than about
///                                 the caller's data: nothing is wrong with the
///                                 waypoints, and a wider Scalar accepts the same
///                                 configuration unchanged. These constructions
///                                 raise it, for the same reason and with the same
///                                 remedy -- carry the command at a scale the type
///                                 spans:
///                                   - a cubic or smoothing spline whose
///                                     polynomial coefficients leave the
///                                     representable range, either by overflowing
///                                     or by flushing a nonzero numerator to zero,
///                                     which drops a term from the polynomial
///                                     rather than rounding it;
///                                   - a smoothing spline whose regularized normal
///                                     equations cannot be squared without
///                                     overflow, which the decomposition that
///                                     solves them does to every entry;
///                                   - a B-spline interpolation whose solved
///                                     control points, or a B-spline whose
///                                     derivative control points, leave the
///                                     representable range.
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
///                                      time. The comparison admits an excess
///                                      of boundary_velocity_rounding_ops * eps
///                                      * v_max, zero by default.
///  * boundary_velocity_rounding_ops_out_of_range
///                                    : boundary_velocity_rounding_ops times eps
///                                      exceeds one. k roundings bound a
///                                      relative error by k u / (1 - k u)
///                                      (@cite higham2002, Lemma 3.1), which the
///                                      admitted excess of k eps covers only
///                                      while k eps <= 1, so a larger count
///                                      names no chain of roundings. Checked
///                                      before either boundary velocity.
///  * unrepresentable_duration        : the arithmetic cannot represent or
///                                      resolve what the answer needs. Five
///                                      cases, all of them facts about the
///                                      arithmetic rather than about the
///                                      kinematics: a duration of the
///                                      constructed profile left the finite
///                                      range; the total underflowed to zero on
///                                      a command with a nonzero displacement,
///                                      which would report an instantaneous
///                                      traversal; a requested retiming differs
///                                      from a duration the profile already
///                                      realizes, or from a coefficient built
///                                      out of that difference, by less than the
///                                      rounding of the expression that would
///                                      solve for it, so it cannot be told apart
///                                      from that duration; the discriminant
///                                      of the shape's own quadratic left the
///                                      finite range, which would drive its root
///                                      to zero and return the shape boundary
///                                      itself; or the constructed profile's
///                                      peak velocity did not come out at or
///                                      below the velocity limit, which, with
///                                      every peak formed at or below it by
///                                      construction, only a non-finite peak
///                                      can produce. Every case
///                                      past the first two says a profile may
///                                      well exist and the arithmetic cannot
///                                      locate it, which is why they are
///                                      reported separately from
///                                      unreachable_duration: the caller should
///                                      change the request, not the limits.
///  * unreachable_boundary_velocity   : the commanded displacement is smaller
///                                      than the distance the fastest admissible
///                                      transition between the two boundary
///                                      velocities already sweeps, so no profile
///                                      of the requested shape realizes it.
///  * acceleration_limit_infeasible   : a trapezoidal command is realizable only
///                                      at an acceleration above the commanded
///                                      limit, under the default
///                                      acceleration_limit_policy::refuse. It
///                                      differs from unreachable_boundary_velocity,
///                                      which no acceleration realizes: here a
///                                      larger limit, or raise_to_feasible, builds it.
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
    boundary_velocity_rounding_ops_out_of_range,
    unrepresentable_duration,
    unreachable_boundary_velocity,
    acceleration_limit_infeasible,
    duration_shorter_than_current,
    unreachable_duration,
};

}

#endif
