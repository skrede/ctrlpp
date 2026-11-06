#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAPEZOIDAL_TYPES_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAPEZOIDAL_TYPES_H

namespace ctrlpp
{

enum class acceleration_limit_policy
{
    refuse,
    raise_to_feasible,
};

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots",
// 2009, Sec. 3.2.7, eq. (3.14), p.71
template <typename Scalar>
struct trapezoidal_disposition
{
    trapezoidal_disposition()
        : commanded_acceleration{}
        , realized_acceleration{}
    {
    }

    Scalar commanded_acceleration;
    Scalar realized_acceleration;
};

}

#endif
