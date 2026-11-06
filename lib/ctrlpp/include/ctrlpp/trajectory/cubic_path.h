#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_PATH_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Table 2.1, p.37; eq. 5.7-5.8

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>

namespace ctrlpp
{

// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
path_point<Scalar> cubic_path(Scalar tau)
{
    auto const tau2 = tau * tau;
    return {
            tau2 * (Scalar{3} - Scalar{2} * tau),
            Scalar{6} * tau * (Scalar{1} - tau),
            Scalar{6} - Scalar{12} * tau,
            Scalar{-12},
    };
}

// @cite biagiotti2009 -- eq. 5.7-5.8
template <typename Scalar>
std::array<Scalar, 3> cubic_path_peak_derivatives()
{
    return {Scalar{1.5}, Scalar{6.0}, Scalar{12.0}};
}

}

#endif
