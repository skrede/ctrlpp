#ifndef HPP_GUARD_CTRLPP_TRAJ_MOTION_LAWS_H
#define HPP_GUARD_CTRLPP_TRAJ_MOTION_LAWS_H

/// @brief Normalized motion law functions for trajectory generation.
///
/// Each law maps normalized time tau in [0,1] to a motion_law_point containing
/// normalized position q in [0,1], velocity dq, acceleration ddq, and jerk dddq.
/// Physical values are obtained via scaling with displacement h and duration T.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Chapter 2

#include "ctrlpp/traj/trajectory_types.h"

#include <cmath>
#include <numbers>

namespace ctrlpp
{

/// @brief Cubic polynomial motion law: q_N(tau) = 3*tau^2 - 2*tau^3
///
/// C1-continuous (zero velocity at endpoints, non-zero acceleration).
///
/// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
auto cubic_law(Scalar tau) -> motion_law_point<Scalar>
{
    auto const tau2 = tau * tau;
    return {
        .q = tau2 * (Scalar{3} - Scalar{2} * tau),
        .dq = Scalar{6} * tau * (Scalar{1} - tau),
        .ddq = Scalar{6} - Scalar{12} * tau,
        .dddq = Scalar{-12},
    };
}

/// @brief Quintic polynomial motion law: q_N(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5
///
/// C2-continuous (zero velocity and acceleration at endpoints).
///
/// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
auto quintic_law(Scalar tau) -> motion_law_point<Scalar>
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

/// @brief Septic polynomial motion law: q_N(tau) = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7
///
/// C3-continuous (zero velocity, acceleration, and jerk at endpoints).
///
/// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
auto septic_law(Scalar tau) -> motion_law_point<Scalar>
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

/// @brief Harmonic motion law: q_N(tau) = (1 - cos(pi*tau)) / 2
///
/// Sinusoidal velocity profile. Non-zero acceleration at endpoints.
///
/// @cite biagiotti2009 -- Sec. 2.2.1, eq. (2.21), p.43
template <typename Scalar>
auto harmonic_law(Scalar tau) -> motion_law_point<Scalar>
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

/// @brief Cycloidal motion law: q_N(tau) = tau - sin(2*pi*tau) / (2*pi)
///
/// Zero acceleration at endpoints. Smooth sinusoidal acceleration profile.
///
/// @cite biagiotti2009 -- Sec. 2.2.2, eq. (2.22), p.44
template <typename Scalar>
auto cycloidal_law(Scalar tau) -> motion_law_point<Scalar>
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

} // namespace ctrlpp

#endif
