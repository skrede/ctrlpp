#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SEPTIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SEPTIC_PATH_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Table 2.1, p.37; eq. 5.11-5.12

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>
#include <cmath>

namespace ctrlpp
{

// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
path_point<Scalar> septic_path(Scalar tau)
{
    auto const tau2 = tau * tau;
    auto const tau3 = tau2 * tau;
    auto const tau4 = tau2 * tau2;
    return {
            tau4 * (Scalar{35} + tau * (Scalar{-84} + tau * (Scalar{70} + tau * Scalar{-20}))),
            tau3 * (Scalar{140} + tau * (Scalar{-420} + tau * (Scalar{420} + tau * Scalar{-140}))),
            tau2 * (Scalar{420} + tau * (Scalar{-1680} + tau * (Scalar{2100} + tau * Scalar{-840}))),
            tau * (Scalar{840} + tau * (Scalar{-5040} + tau * (Scalar{8400} + tau * Scalar{-4200}))),
    };
}

// ddq_max ~ 7.5132 (non-trivial closed form involving nested radicals).
// @cite biagiotti2009 -- eq. 5.11-5.12
template <typename Scalar>
std::array<Scalar, 3> septic_path_peak_derivatives()
{
    return {
            Scalar{35} / Scalar{16},
            static_cast<Scalar>(7.5132),
            Scalar{52.5},
    };
}

}

#endif
