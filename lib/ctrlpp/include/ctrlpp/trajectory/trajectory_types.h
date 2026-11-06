#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_TYPES_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_TYPES_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 5.2.1, eq. (5.5)-(5.6), p.230

#include "ctrlpp/trajectory/trajectory_error.h"

#include "ctrlpp/types.h"

#include <cstddef>

namespace ctrlpp
{

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 1.6, p.10
// An aggregate without a constructor: designated initialization stores each field as it is computed, where a
// member-wise constructor holds every argument until the call and grows the per-cycle evaluate frames.
template <typename Scalar, std::size_t ND>
struct trajectory_point
{
    Vector<Scalar, ND> position;
    Vector<Scalar, ND> velocity;
    Vector<Scalar, ND> acceleration;
};

template <typename Scalar>
struct path_point
{
    path_point()
        : q{}
        , dq{}
        , ddq{}
        , dddq{}
    {
    }

    path_point(Scalar q_value, Scalar dq_value, Scalar ddq_value, Scalar dddq_value)
        : q{q_value}
        , dq{dq_value}
        , ddq{ddq_value}
        , dddq{dddq_value}
    {
    }

    Scalar q;
    Scalar dq;
    Scalar ddq;
    Scalar dddq;
};

}

#endif
