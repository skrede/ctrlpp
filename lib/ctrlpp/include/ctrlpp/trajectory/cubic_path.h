#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_PATH_H

/// @brief Cubic polynomial path: q_N(tau) = 3*tau^2 - 2*tau^3
///
/// C1-continuous (zero velocity at endpoints, non-zero acceleration).
/// Peak normalized derivatives: dq_max = 1.5, ddq_max = 6.0, dddq_max = 12.0.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Table 2.1, p.37; eq. 5.7-5.8

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>

namespace ctrlpp
{

/// @brief Cubic polynomial path: q_N(tau) = 3*tau^2 - 2*tau^3
///
/// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
auto cubic_path(Scalar tau) -> path_point<Scalar>
{
    auto const tau2 = tau * tau;
    return {
        .q = tau2 * (Scalar{3} - Scalar{2} * tau),
        .dq = Scalar{6} * tau * (Scalar{1} - tau),
        .ddq = Scalar{6} - Scalar{12} * tau,
        .dddq = Scalar{-12},
    };
}

/// @brief Peak normalized derivatives for cubic path: {dq_max, ddq_max, dddq_max}.
///
/// @cite biagiotti2009 -- eq. 5.7-5.8
template <typename Scalar>
auto cubic_path_peak_derivatives() -> std::array<Scalar, 3>
{
    return {Scalar{1.5}, Scalar{6.0}, Scalar{12.0}};
}

} // namespace ctrlpp

#endif
