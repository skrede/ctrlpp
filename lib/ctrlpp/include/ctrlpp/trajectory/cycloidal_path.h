#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_CYCLOIDAL_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_CYCLOIDAL_PATH_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Sec. 2.2.2, eq. (2.22), p.43-44; p.235

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>
#include <cmath>
#include <numbers>

namespace ctrlpp
{

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 2.2.2, eq. (2.22), p.43-44
template <typename Scalar>
path_point<Scalar> cycloidal_path(Scalar tau)
{
    auto const pi = std::numbers::pi_v<Scalar>;
    auto const two_pi = Scalar{2} * pi;
    auto const two_pi_tau = two_pi * tau;
    auto const sin_val = std::sin(two_pi_tau);
    auto const cos_val = std::cos(two_pi_tau);
    auto const pi2 = pi * pi;
    return {
            tau - sin_val / two_pi,
            Scalar{1} - cos_val,
            two_pi * sin_val,
            Scalar{4} * pi2 * cos_val,
    };
}

// @cite biagiotti2009 -- p.235
template <typename Scalar>
std::array<Scalar, 3> cycloidal_path_peak_derivatives()
{
    auto const pi = std::numbers::pi_v<Scalar>;
    auto const pi2 = pi * pi;
    return {Scalar{2}, Scalar{2} * pi, Scalar{4} * pi2};
}

}

#endif
