#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_HARMONIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_HARMONIC_PATH_H

/// @brief Harmonic path: q_N(tau) = (1 - cos(pi*tau)) / 2
///
/// Sinusoidal velocity profile. Non-zero acceleration at endpoints.
/// Peak normalized derivatives: dq_max = pi/2, ddq_max = pi^2/2, dddq_max = pi^3/2.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 2.2.1, eq. (2.21), p.43; p.236

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>
#include <cmath>
#include <numbers>

namespace ctrlpp
{

/// @brief Harmonic path: q_N(tau) = (1 - cos(pi*tau)) / 2
///
/// @cite biagiotti2009 -- Sec. 2.2.1, eq. (2.21), p.43
template <typename Scalar>
auto harmonic_path(Scalar tau) -> path_point<Scalar>
{
    auto const pi = std::numbers::pi_v<Scalar>;
    auto const pi_tau = pi * tau;
    auto const sin_val = std::sin(pi_tau);
    auto const cos_val = std::cos(pi_tau);
    auto const pi2 = pi * pi;
    auto const pi3 = pi2 * pi;
    return {
        .q = (Scalar{1} - cos_val) / Scalar{2},
        .dq = pi / Scalar{2} * sin_val,
        .ddq = pi2 / Scalar{2} * cos_val,
        .dddq = -pi3 / Scalar{2} * sin_val,
    };
}

/// @brief Peak normalized derivatives for harmonic path: {dq_max, ddq_max, dddq_max}.
///
/// @cite biagiotti2009 -- p.236
template <typename Scalar>
auto harmonic_path_peak_derivatives() -> std::array<Scalar, 3>
{
    auto const pi = std::numbers::pi_v<Scalar>;
    auto const pi2 = pi * pi;
    auto const pi3 = pi2 * pi;
    return {pi / Scalar{2}, pi2 / Scalar{2}, pi3 / Scalar{2}};
}

} // namespace ctrlpp

#endif
