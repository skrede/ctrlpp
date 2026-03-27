#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_TRAJECTORY_H

/// @brief Cubic polynomial trajectory segment with arbitrary velocity BCs.
///
/// Interpolates between two points with specified endpoint velocities using a
/// degree-3 polynomial in normalized time tau = t/T. Coefficients stored in
/// normalized basis to prevent ill-conditioning.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 2.1.4, eq. (2.2), p.24

#include "ctrlpp/trajectory/detail/polynomial_eval.h"
#include "ctrlpp/trajectory/trajectory_segment.h"
#include "ctrlpp/trajectory/trajectory_types.h"

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctrlpp
{

/// @brief Cubic polynomial trajectory segment: q(tau) = c0 + c1*tau + c2*tau^2 + c3*tau^3
///
/// @cite biagiotti2009 -- Sec. 2.1.4, eq. (2.2), p.24: four-coefficient cubic in normalized time
template <typename Scalar, std::size_t ND>
struct cubic_trajectory
{
    std::array<Vector<Scalar, ND>, 4> coeffs; ///< Normalized polynomial coefficients
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

/// @brief Create cubic trajectory from boundary conditions q(0)=q0, q(T)=q1, dq(0)=v0, dq(T)=v1.
///
/// @cite biagiotti2009 -- Sec. 2.1.4, eq. (2.2), p.24: coefficient derivation from BCs
/// Coefficients derived from B&M eq. (2.2) in normalized time tau = t/T.
template <typename Scalar, int Rows>
auto make_cubic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,
    Eigen::Matrix<Scalar, Rows, 1> const& q1,
    Eigen::Matrix<Scalar, Rows, 1> const& v0,
    Eigen::Matrix<Scalar, Rows, 1> const& v1,
    Scalar duration) -> cubic_trajectory<Scalar, static_cast<std::size_t>(Rows)>
{
    constexpr auto ND = static_cast<std::size_t>(Rows);
    auto const T = duration;
    Vector<Scalar, ND> const h = (q1 - q0).eval();
    Vector<Scalar, ND> const c0 = q0;
    Vector<Scalar, ND> const c1 = (v0 * T).eval();
    // @cite biagiotti2009 -- Sec. 2.1.4, eq. (2.2), p.24: c2, c3 from BC linear system
    Vector<Scalar, ND> const c2 = (Scalar{3} * h - T * (Scalar{2} * v0 + v1)).eval();
    Vector<Scalar, ND> const c3 = (Scalar{-2} * h + T * (v0 + v1)).eval();
    return {.coeffs = {c0, c1, c2, c3}, .dur = duration};
}

static_assert(trajectory_segment<cubic_trajectory<double, 1>, double, 1>);

} // namespace ctrlpp

#endif
