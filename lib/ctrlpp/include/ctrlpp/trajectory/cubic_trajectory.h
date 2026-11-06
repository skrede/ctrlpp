#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_TRAJECTORY_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Sec. 2.1.4, eq. (2.2), p.24

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/trajectory/detail/polynomial_eval.h"

#include "ctrlpp/expected.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 2.1.4, eq. (2.2), p.24: four-coefficient cubic in normalized time
template <typename Scalar, std::size_t ND>
struct cubic_trajectory
{
    std::array<Vector<Scalar, ND>, 4> coeffs; // in normalized time tau = t / dur
    Scalar dur;                               // [s]

    trajectory_point<Scalar, ND> evaluate(Scalar t) const
    {
        auto const tau = std::clamp(t, Scalar{0}, dur) / dur;
        auto const inv_T = Scalar{1} / dur;
        return {
                .position = detail::horner_eval<Scalar, ND>(coeffs, tau),
                .velocity = detail::horner_deriv1<Scalar, ND>(coeffs, tau) * inv_T,
                .acceleration = detail::horner_deriv2<Scalar, ND>(coeffs, tau) * (inv_T * inv_T),
        };
    }

    Scalar duration() const { return dur; }
};

// @cite biagiotti2009 -- Sec. 2.1.4, eq. (2.2), p.24: coefficient derivation from BCs
template <typename Scalar, int Rows>
ctrlpp::expected<cubic_trajectory<Scalar, static_cast<std::size_t>(Rows)>, trajectory_error> make_cubic_trajectory(Eigen::Matrix<Scalar, Rows, 1> const& q0,
        Eigen::Matrix<Scalar, Rows, 1> const& q1, Eigen::Matrix<Scalar, Rows, 1> const& v0, Eigen::Matrix<Scalar, Rows, 1> const& v1, Scalar duration)
{
    constexpr auto ND = static_cast<std::size_t>(Rows);
    if(!std::isfinite(duration) || !q0.allFinite() || !q1.allFinite() || !v0.allFinite() || !v1.allFinite())
    {
        return ctrlpp::unexpected(trajectory_error::non_finite_input);
    }
    if(duration <= Scalar{0})
    {
        return ctrlpp::unexpected(trajectory_error::non_positive_duration);
    }
    auto const T = duration;
    Vector<Scalar, ND> const h = (q1 - q0).eval();
    Vector<Scalar, ND> const c0 = q0;
    Vector<Scalar, ND> const c1 = (v0 * T).eval();
    // @cite biagiotti2009 -- Sec. 2.1.4, eq. (2.2), p.24: c2, c3 from BC linear system
    Vector<Scalar, ND> const c2 = (Scalar{3} * h - T * (Scalar{2} * v0 + v1)).eval();
    Vector<Scalar, ND> const c3 = (Scalar{-2} * h + T * (v0 + v1)).eval();
    return cubic_trajectory<Scalar, ND>{.coeffs = {c0, c1, c2, c3}, .dur = duration};
}

static_assert(trajectory_segment<cubic_trajectory<double, 1>, double, 1>);

}

#endif
