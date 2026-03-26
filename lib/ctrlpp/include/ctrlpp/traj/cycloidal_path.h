#ifndef HPP_GUARD_CTRLPP_TRAJ_CYCLOIDAL_PATH_H
#define HPP_GUARD_CTRLPP_TRAJ_CYCLOIDAL_PATH_H

/// @brief Cycloidal path: q_N(tau) = tau - sin(2*pi*tau) / (2*pi)
///
/// Zero acceleration at endpoints. Smooth sinusoidal acceleration profile.
/// Peak normalized derivatives: dq_max = 2.0, ddq_max = 2*pi, dddq_max = 4*pi^2.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 2.2.2, eq. (2.22), p.44; p.235

#include "ctrlpp/traj/trajectory_types.h"

#include <array>
#include <cmath>
#include <numbers>

namespace ctrlpp
{

/// @brief Cycloidal path: q_N(tau) = tau - sin(2*pi*tau) / (2*pi)
///
/// @cite biagiotti2009 -- Sec. 2.2.2, eq. (2.22), p.44
template <typename Scalar>
auto cycloidal_path(Scalar tau) -> path_point<Scalar>
{
    auto const pi = std::numbers::pi_v<Scalar>;
    auto const two_pi = Scalar{2} * pi;
    auto const two_pi_tau = two_pi * tau;
    auto const sin_val = std::sin(two_pi_tau);
    auto const cos_val = std::cos(two_pi_tau);
    auto const pi2 = pi * pi;
    return {
        .q = tau - sin_val / two_pi,
        .dq = Scalar{1} - cos_val,
        .ddq = two_pi * sin_val,
        .dddq = Scalar{4} * pi2 * cos_val,
    };
}

/// @brief Peak normalized derivatives for cycloidal path: {dq_max, ddq_max, dddq_max}.
///
/// @cite biagiotti2009 -- p.235
template <typename Scalar>
auto cycloidal_path_peak_derivatives() -> std::array<Scalar, 3>
{
    auto const pi = std::numbers::pi_v<Scalar>;
    auto const pi2 = pi * pi;
    return {Scalar{2}, Scalar{2} * pi, Scalar{4} * pi2};
}

} // namespace ctrlpp

#endif
