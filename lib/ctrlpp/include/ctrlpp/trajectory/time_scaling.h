#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TIME_SCALING_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TIME_SCALING_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 5.2.1, eq. (5.5)-(5.6), p.230-231

#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 5.2.1, eq. (5.6), p.230
template <typename Scalar>
Scalar compute_min_duration(Scalar h, std::array<Scalar, 3> const& peak_derivs, Scalar v_max, Scalar a_max, Scalar j_max)
{
    auto const abs_h = std::abs(h);
    auto const T_vel = abs_h * peak_derivs[0] / v_max;
    auto const T_acc = std::sqrt(abs_h * peak_derivs[1] / a_max);
    auto const T_jrk = std::cbrt(abs_h * peak_derivs[2] / j_max);
    return std::max({T_vel, T_acc, T_jrk});
}

// @cite biagiotti2009 -- Sec. 5.2.1, eq. (5.6), p.230
// @cite biagiotti2009 -- Sec. 5.2.2 -- per-DOF timing extends scalar formulation to vector case
template <typename Scalar, std::size_t N>
std::array<Scalar, N> compute_min_duration(std::array<Scalar, N> const& h, std::array<Scalar, 3> const& peak_derivs, Scalar v_max, Scalar a_max, Scalar j_max)
{
    std::array<Scalar, N> result{};
    for(std::size_t i = 0; i < N; ++i)
    {
        result[i] = compute_min_duration(h[i], peak_derivs, v_max, a_max, j_max);
    }
    return result;
}

// @cite biagiotti2009 -- Sec. 5.2.1, eq. (5.5)-(5.6), p.230 -- kinematic limit satisfaction
template <typename Scalar, std::size_t N>
Scalar compute_min_duration_sync(std::array<Scalar, N> const& h, std::array<Scalar, 3> const& peak_derivs, Scalar v_max, Scalar a_max, Scalar j_max)
{
    auto const per_dof = compute_min_duration(h, peak_derivs, v_max, a_max, j_max);
    return *std::max_element(per_dof.begin(), per_dof.end());
}

}

#endif
