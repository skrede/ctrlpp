#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_SEGMENT_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_SEGMENT_H

/// @brief Concept constraining trajectory segment types.
///
/// A trajectory segment must provide evaluate(t) returning a trajectory_point
/// and duration() returning a scalar. Uses named method pattern (not operator()).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009

#include "ctrlpp/trajectory/trajectory_types.h"

#include <concepts>
#include <cstddef>

namespace ctrlpp
{

/// @brief Constrains a type to be a trajectory segment with evaluate(t) and duration().
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Sec. 2.1.1 -- segment evaluation interface q(t) -> (pos, vel, acc)
template <typename S, typename Scalar, std::size_t ND>
concept trajectory_segment = requires(const S& seg, Scalar t) {
    { seg.evaluate(t) } -> std::convertible_to<trajectory_point<Scalar, ND>>;
    { seg.duration() } -> std::convertible_to<Scalar>;
};

} // namespace ctrlpp

#endif
