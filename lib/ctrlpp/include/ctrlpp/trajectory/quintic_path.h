#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_QUINTIC_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_QUINTIC_PATH_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Table 2.1, p.37; eq. 5.9-5.10

#include "ctrlpp/trajectory/trajectory_types.h"

#include <array>
#include <cmath>

namespace ctrlpp
{

// @cite biagiotti2009 -- Table 2.1, p.37
template <typename Scalar>
path_point<Scalar> quintic_path(Scalar tau)
{
    auto const tau2 = tau * tau;
    auto const tau3 = tau2 * tau;
    return {
            tau3 * (Scalar{10} + tau * (Scalar{-15} + tau * Scalar{6})),
            tau2 * (Scalar{30} + tau * (Scalar{-60} + tau * Scalar{30})),
            tau * (Scalar{60} + tau * (Scalar{-180} + tau * Scalar{120})),
            Scalar{60} + tau * (Scalar{-360} + tau * Scalar{360}),
    };
}

// @cite biagiotti2009 -- eq. 5.9-5.10
template <typename Scalar>
std::array<Scalar, 3> quintic_path_peak_derivatives()
{
    using std::sqrt;
    return {
            Scalar{15} / Scalar{8},
            Scalar{10} * static_cast<Scalar>(sqrt(Scalar{3})) / Scalar{3},
            Scalar{60},
    };
}

}

#endif
