#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_QUINTIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_QUINTIC_PATH_H

/// @brief Quintic polynomial path: q_N(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5
///
/// C2-continuous (zero velocity and acceleration at endpoints).
/// Peak normalized derivatives: dq_max = 15/8, ddq_max = 10*sqrt(3)/3, dddq_max = 60.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Table 2.1, p.37; eq. 5.9-5.10

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>
#include <cmath>

namespace ctrlpp
{

/// @brief Quintic polynomial path: q_N(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5
///
/// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
auto quintic_path(Scalar tau) -> path_point<Scalar>
{
    auto const tau2 = tau * tau;
    auto const tau3 = tau2 * tau;
    return {
        .q = tau3 * (Scalar{10} + tau * (Scalar{-15} + tau * Scalar{6})),
        .dq = tau2 * (Scalar{30} + tau * (Scalar{-60} + tau * Scalar{30})),
        .ddq = tau * (Scalar{60} + tau * (Scalar{-180} + tau * Scalar{120})),
        .dddq = Scalar{60} + tau * (Scalar{-360} + tau * Scalar{360}),
    };
}

/// @brief Peak normalized derivatives for quintic path: {dq_max, ddq_max, dddq_max}.
///
/// ddq_max = 10*sqrt(3)/3 is not constexpr due to sqrt.
///
/// @cite biagiotti2009 -- eq. 5.9-5.10
template <typename Scalar>
auto quintic_path_peak_derivatives() -> std::array<Scalar, 3>
{
    using std::sqrt;
    return {
        Scalar{15} / Scalar{8},
        Scalar{10} * static_cast<Scalar>(sqrt(Scalar{3})) / Scalar{3},
        Scalar{60},
    };
}

} // namespace ctrlpp

#endif
