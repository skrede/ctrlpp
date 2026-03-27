#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_QUINTIC_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_QUINTIC_TRAJECTORY_H

/// @brief Quintic polynomial trajectory segment with arbitrary velocity and acceleration BCs.
///
/// Interpolates between two points with specified endpoint velocities and
/// accelerations using a degree-5 polynomial in normalized time tau = t/T.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 2.1.5, eq. (2.5), p.27

#include "ctrlpp/trajectory/detail/polynomial_eval.h"
#include "ctrlpp/trajectory/trajectory_segment.h"
#include "ctrlpp/trajectory/trajectory_types.h"

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctrlpp
{

/// @brief Quintic polynomial trajectory segment: c0 + c1*tau + ... + c5*tau^5
///
/// @cite biagiotti2009 -- Sec. 2.1.5, eq. (2.5), p.27: six-coefficient quintic in normalized time
template <typename Scalar, std::size_t ND>
struct quintic_trajectory
{
    std::array<Vector<Scalar, ND>, 6> coeffs; ///< Normalized polynomial coefficients
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

/// @brief Create quintic trajectory from boundary conditions q, dq, ddq at both endpoints.
///
/// @cite biagiotti2009 -- Sec. 2.1.5, eq. (2.5), p.27: coefficient derivation from 6 BCs
/// Coefficients derived from B&M eq. (2.5) in normalized time tau = t/T.
template <typename Scalar, int Rows>
auto make_quintic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,
    Eigen::Matrix<Scalar, Rows, 1> const& q1,
    Eigen::Matrix<Scalar, Rows, 1> const& v0,
    Eigen::Matrix<Scalar, Rows, 1> const& v1,
    Eigen::Matrix<Scalar, Rows, 1> const& a0,
    Eigen::Matrix<Scalar, Rows, 1> const& a1,
    Scalar duration) -> quintic_trajectory<Scalar, static_cast<std::size_t>(Rows)>
{
    constexpr auto ND = static_cast<std::size_t>(Rows);
    auto const T = duration;
    auto const T2 = T * T;
    Vector<Scalar, ND> const h = (q1 - q0).eval();

    Vector<Scalar, ND> const c0 = q0;
    Vector<Scalar, ND> const c1 = (v0 * T).eval();
    Vector<Scalar, ND> const c2 = (a0 * T2 / Scalar{2}).eval();
    // @cite biagiotti2009 -- Sec. 2.1.5, eq. (2.5), p.27: c3-c5 from BC linear system
    Vector<Scalar, ND> const c3 = ((Scalar{20} * h - (Scalar{8} * v1 + Scalar{12} * v0) * T
                                     - (Scalar{3} * a0 - a1) * T2)
                                    / Scalar{2})
                                      .eval();
    Vector<Scalar, ND> const c4 = ((Scalar{-30} * h + (Scalar{14} * v1 + Scalar{16} * v0) * T
                                     + (Scalar{3} * a0 - Scalar{2} * a1) * T2)
                                    / Scalar{2})
                                      .eval();
    Vector<Scalar, ND> const c5 =
        ((Scalar{12} * h - Scalar{6} * (v1 + v0) * T + (a1 - a0) * T2) / Scalar{2}).eval();

    return {.coeffs = {c0, c1, c2, c3, c4, c5}, .dur = duration};
}

static_assert(trajectory_segment<quintic_trajectory<double, 1>, double, 1>);

} // namespace ctrlpp

#endif
