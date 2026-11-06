#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_PATH_SEGMENT_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_PATH_SEGMENT_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Ch. 2 (motion laws)

#include "ctrlpp/trajectory/trajectory_types.h"

#include <concepts>

namespace ctrlpp
{

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
//   Machines and Robots", 2008, Sec. 5.2.1 -- motion law abstraction tau in [0,1]
template <typename S, typename Scalar>
concept path_segment = requires(const S& seg, Scalar tau) {
    { seg.evaluate(tau) } -> std::convertible_to<path_point<Scalar>>;
    { seg.duration() } -> std::convertible_to<Scalar>;
};

}

#endif
