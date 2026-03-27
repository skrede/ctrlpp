#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SEPTIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SEPTIC_PATH_H

/// @brief Septic polynomial path: q_N(tau) = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7
///
/// C3-continuous (zero velocity, acceleration, and jerk at endpoints).
/// Peak normalized derivatives: dq_max = 35/16, ddq_max ~ 7.5132, dddq_max = 52.5.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Table 2.1, p.37; eq. 5.11-5.12

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>
#include <cmath>

namespace ctrlpp
{

/// @brief Septic polynomial path: q_N(tau) = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7
///
/// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
auto septic_path(Scalar tau) -> path_point<Scalar>
{
    auto const tau2 = tau * tau;
    auto const tau3 = tau2 * tau;
    auto const tau4 = tau2 * tau2;
    return {
        .q = tau4 * (Scalar{35} + tau * (Scalar{-84} + tau * (Scalar{70} + tau * Scalar{-20}))),
        .dq = tau3 * (Scalar{140} + tau * (Scalar{-420} + tau * (Scalar{420} + tau * Scalar{-140}))),
        .ddq = tau2 * (Scalar{420} + tau * (Scalar{-1680} + tau * (Scalar{2100} + tau * Scalar{-840}))),
        .dddq = tau * (Scalar{840} + tau * (Scalar{-5040} + tau * (Scalar{8400} + tau * Scalar{-4200}))),
    };
}

/// @brief Peak normalized derivatives for septic path: {dq_max, ddq_max, dddq_max}.
///
/// ddq_max ~ 7.5132 (non-trivial closed form involving nested radicals).
///
/// @cite biagiotti2009 -- eq. 5.11-5.12
template <typename Scalar>
auto septic_path_peak_derivatives() -> std::array<Scalar, 3>
{
    return {
        Scalar{35} / Scalar{16},
        static_cast<Scalar>(7.5132),
        Scalar{52.5},
    };
}

} // namespace ctrlpp

#endif
