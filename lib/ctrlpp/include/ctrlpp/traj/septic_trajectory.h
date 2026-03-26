#ifndef HPP_GUARD_CTRLPP_TRAJ_SEPTIC_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJ_SEPTIC_TRAJECTORY_H

/// @brief Septic polynomial trajectory segment with arbitrary BCs up to jerk.
///
/// Interpolates between two points with specified endpoint velocities,
/// accelerations, and jerks using a degree-7 polynomial in normalized time tau = t/T.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 2.1.6, eq. (2.6), p.29

#include "ctrlpp/traj/detail/polynomial_eval.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctrlpp
{

/// @brief Septic polynomial trajectory segment: c0 + c1*tau + ... + c7*tau^7
template <typename Scalar, std::size_t ND>
struct septic_trajectory
{
    std::array<Vector<Scalar, ND>, 8> coeffs; ///< Normalized polynomial coefficients
    Scalar dur;                                ///< Segment duration [s]

    auto evaluate(Scalar t) const -> trajectory_point<Scalar, ND>
    {
        auto const tau = std::clamp(t, Scalar{0}, dur) / dur;
        auto const inv_T = Scalar{1} / dur;
        return {
            .position = detail::horner_eval<Scalar, ND>(coeffs, tau),
            .velocity = detail::horner_deriv1<Scalar, ND>(coeffs, tau) * inv_T,
            .acceleration = detail::horner_deriv2<Scalar, ND>(coeffs, tau) * (inv_T * inv_T),
        };
    }

    auto duration() const -> Scalar { return dur; }
};

/// @brief Create septic trajectory from boundary conditions q, dq, ddq, d3q at both endpoints.
///
/// Coefficients derived from B&M eq. (2.6) in normalized time tau = t/T.
template <typename Scalar, int Rows>
auto make_septic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,
    Eigen::Matrix<Scalar, Rows, 1> const& q1,
    Eigen::Matrix<Scalar, Rows, 1> const& v0,
    Eigen::Matrix<Scalar, Rows, 1> const& v1,
    Eigen::Matrix<Scalar, Rows, 1> const& a0,
    Eigen::Matrix<Scalar, Rows, 1> const& a1,
    Eigen::Matrix<Scalar, Rows, 1> const& j0,
    Eigen::Matrix<Scalar, Rows, 1> const& j1,
    Scalar duration) -> septic_trajectory<Scalar, static_cast<std::size_t>(Rows)>
{
    constexpr auto ND = static_cast<std::size_t>(Rows);
    auto const T = duration;
    auto const T2 = T * T;
    auto const T3 = T2 * T;
    Vector<Scalar, ND> const h = (q1 - q0).eval();

    Vector<Scalar, ND> const c0 = q0;
    Vector<Scalar, ND> const c1 = (v0 * T).eval();
    Vector<Scalar, ND> const c2 = (a0 * T2 / Scalar{2}).eval();
    Vector<Scalar, ND> const c3 = (j0 * T3 / Scalar{6}).eval();
    Vector<Scalar, ND> const c4 =
        ((Scalar{210} * h
          - T * ((Scalar{30} * a0 - Scalar{15} * a1) * T + (Scalar{4} * j0 + j1) * T2
                 + Scalar{120} * v0 + Scalar{90} * v1))
         / Scalar{6})
            .eval();
    Vector<Scalar, ND> const c5 =
        ((Scalar{-168} * h
          + T * ((Scalar{20} * a0 - Scalar{14} * a1) * T + (Scalar{2} * j0 + j1) * T2
                 + Scalar{90} * v0 + Scalar{78} * v1))
         / Scalar{2})
            .eval();
    Vector<Scalar, ND> const c6 =
        ((Scalar{420} * h
          - T * ((Scalar{45} * a0 - Scalar{39} * a1) * T
                 + (Scalar{4} * j0 + Scalar{3} * j1) * T2 + Scalar{216} * v0
                 + Scalar{204} * v1))
         / Scalar{6})
            .eval();
    Vector<Scalar, ND> const c7 =
        ((Scalar{-120} * h
          + T * ((Scalar{12} * a0 - Scalar{12} * a1) * T + (j0 + j1) * T2 + Scalar{60} * v0
                 + Scalar{60} * v1))
         / Scalar{6})
            .eval();

    return {.coeffs = {c0, c1, c2, c3, c4, c5, c6, c7}, .dur = duration};
}

static_assert(trajectory_segment<septic_trajectory<double, 1>, double, 1>);

} // namespace ctrlpp

#endif
