#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_SEGMENT_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_SEGMENT_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2009

#include "ctrlpp/trajectory/trajectory_types.h"

#include <cstddef>
#include <concepts>

namespace ctrlpp
{

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
//   Machines and Robots", 2009, Sec. 1.6, p.10
template <typename S, typename Scalar, std::size_t ND>
concept trajectory_segment = requires(const S& seg, Scalar t) {
    { seg.evaluate(t) } -> std::convertible_to<trajectory_point<Scalar, ND>>;
    { seg.duration() } -> std::convertible_to<Scalar>;
};

}

#endif
