#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SEPTIC_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SEPTIC_TRAJECTORY_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Sec. 2.1.6, eq. (2.6), p.28-29

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

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 2.1.6, eq. (2.6), p.28
template <typename Scalar, std::size_t ND>
struct septic_trajectory
{
    std::array<Vector<Scalar, ND>, 8> coeffs; // in normalized time tau = t / dur
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

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 2.1.6, p.29
template <typename Scalar, int Rows>
ctrlpp::expected<septic_trajectory<Scalar, static_cast<std::size_t>(Rows)>, trajectory_error> make_septic_trajectory(Eigen::Matrix<Scalar, Rows, 1> const& q0,
        Eigen::Matrix<Scalar, Rows, 1> const& q1, Eigen::Matrix<Scalar, Rows, 1> const& v0, Eigen::Matrix<Scalar, Rows, 1> const& v1, Eigen::Matrix<Scalar, Rows, 1> const& a0,
        Eigen::Matrix<Scalar, Rows, 1> const& a1, Eigen::Matrix<Scalar, Rows, 1> const& j0, Eigen::Matrix<Scalar, Rows, 1> const& j1, Scalar duration)
{
    if(!std::isfinite(duration) || !q0.allFinite() || !q1.allFinite() || !v0.allFinite() || !v1.allFinite() || !a0.allFinite() || !a1.allFinite() || !j0.allFinite() || !j1.allFinite())
    {
        return ctrlpp::unexpected(trajectory_error::non_finite_input);
    }
    if(duration <= Scalar{0})
    {
        return ctrlpp::unexpected(trajectory_error::non_positive_duration);
    }
    auto const T = duration;
    auto const T2 = T * T;
    auto const T3 = T2 * T;
    auto const h = (q1 - q0).eval();

    auto const c0 = q0;
    auto const c1 = (v0 * T).eval();
    auto const c2 = (a0 * T2 / Scalar{2}).eval();
    auto const c3 = (j0 * T3 / Scalar{6}).eval();
    auto const c4 = ((Scalar{210} * h - T * ((Scalar{30} * a0 - Scalar{15} * a1) * T + (Scalar{4} * j0 + j1) * T2 + Scalar{120} * v0 + Scalar{90} * v1)) / Scalar{6}).eval();
    auto const c5 = ((Scalar{-168} * h + T * ((Scalar{20} * a0 - Scalar{14} * a1) * T + (Scalar{2} * j0 + j1) * T2 + Scalar{90} * v0 + Scalar{78} * v1)) / Scalar{2}).eval();
    auto const c6 =
            ((Scalar{420} * h - T * ((Scalar{45} * a0 - Scalar{39} * a1) * T + (Scalar{4} * j0 + Scalar{3} * j1) * T2 + Scalar{216} * v0 + Scalar{204} * v1)) / Scalar{6}).eval();
    auto const c7 = ((Scalar{-120} * h + T * ((Scalar{12} * a0 - Scalar{12} * a1) * T + (j0 + j1) * T2 + Scalar{60} * v0 + Scalar{60} * v1)) / Scalar{6}).eval();

    return septic_trajectory<Scalar, static_cast<std::size_t>(Rows)>{.coeffs = {c0, c1, c2, c3, c4, c5, c6, c7}, .dur = duration};
}

static_assert(trajectory_segment<septic_trajectory<double, 1>, double, 1>);

}

#endif
