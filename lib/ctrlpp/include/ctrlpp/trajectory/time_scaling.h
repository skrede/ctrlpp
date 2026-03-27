#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TIME_SCALING_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TIME_SCALING_H

/// @brief Kinematic time scaling for normalized paths.
///
/// Computes minimum trajectory duration T from displacement magnitude h,
/// peak normalized path derivatives, and kinematic limits (v_max, a_max, j_max).
/// T = max(h*dq_max/v_max, sqrt(h*ddq_max/a_max), cbrt(h*dddq_max/j_max)).
///
/// Provides scalar (single DOF), vector (per-DOF), and synchronized (max of
/// per-DOF) overloads. For use with elementary paths only -- composite profiles
/// (trapezoidal, double-S) solve timing internally (D-14).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 5.2.1, eq. (5.5)-(5.6), p.230-231

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace ctrlpp
{

/// @brief Compute minimum duration T for a single DOF.
///
/// @param h Displacement magnitude (|q1 - q0|).
/// @param peak_derivs Peak normalized derivatives {dq_max, ddq_max, dddq_max}
///        from the path's *_peak_derivatives() function.
/// @param v_max Maximum velocity limit.
/// @param a_max Maximum acceleration limit.
/// @param j_max Maximum jerk limit.
/// @return Minimum feasible duration T [s].
///
/// @cite biagiotti2009 -- Sec. 5.2.1, eq. (5.6), p.230
template <typename Scalar>
auto compute_min_duration(
    Scalar h,
    std::array<Scalar, 3> const& peak_derivs,
    Scalar v_max,
    Scalar a_max,
    Scalar j_max) -> Scalar
{
    auto const abs_h = std::abs(h);
    auto const T_vel = abs_h * peak_derivs[0] / v_max;
    auto const T_acc = std::sqrt(abs_h * peak_derivs[1] / a_max);
    auto const T_jrk = std::cbrt(abs_h * peak_derivs[2] / j_max);
    return std::max({T_vel, T_acc, T_jrk});
}

/// @brief Compute per-DOF minimum durations for N degrees of freedom.
///
/// Returns an array of N durations, one per DOF. Each is independently
/// computed from the same path type and kinematic limits.
///
/// @cite biagiotti2009 -- Sec. 5.2.1, eq. (5.6), p.230
/// @cite biagiotti2009 -- Sec. 5.2.2 -- per-DOF timing extends scalar formulation to vector case
template <typename Scalar, std::size_t N>
auto compute_min_duration(
    std::array<Scalar, N> const& h,
    std::array<Scalar, 3> const& peak_derivs,
    Scalar v_max,
    Scalar a_max,
    Scalar j_max) -> std::array<Scalar, N>
{
    std::array<Scalar, N> result{};
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = compute_min_duration(h[i], peak_derivs, v_max, a_max, j_max);
    }
    return result;
}

/// @brief Compute synchronized (maximum) minimum duration across all DOFs.
///
/// Returns a single T suitable for synchronizing all joints to complete
/// their motion simultaneously. Minimum-time computation satisfying kinematic limits.
///
/// @cite biagiotti2009 -- Sec. 5.2.1, eq. (5.5)-(5.6), p.230 -- kinematic limit satisfaction
template <typename Scalar, std::size_t N>
auto compute_min_duration_sync(
    std::array<Scalar, N> const& h,
    std::array<Scalar, 3> const& peak_derivs,
    Scalar v_max,
    Scalar a_max,
    Scalar j_max) -> Scalar
{
    auto const per_dof = compute_min_duration(h, peak_derivs, v_max, a_max, j_max);
    return *std::max_element(per_dof.begin(), per_dof.end());
}

} // namespace ctrlpp

#endif
