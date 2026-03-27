#ifndef HPP_GUARD_CTRLPP_TRAJ_PATH_SEGMENT_H
#define HPP_GUARD_CTRLPP_TRAJ_PATH_SEGMENT_H

/// @brief Concept constraining normalized path segment types.
///
/// A path segment must provide evaluate(tau) returning a path_point
/// and duration() returning a scalar. Paths operate in normalized time
/// tau in [0,1] and produce scalar path_point outputs.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Ch. 2 (motion laws)

#include "ctrlpp/traj/trajectory_types.h"

#include <concepts>

namespace ctrlpp
{

/// @brief Constrains a type to be a path segment with evaluate(tau) and duration().
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Sec. 5.2.1 -- motion law abstraction tau in [0,1]
template <typename S, typename Scalar>
concept path_segment = requires(const S& seg, Scalar tau) {
    { seg.evaluate(tau) } -> std::convertible_to<path_point<Scalar>>;
    { seg.duration() } -> std::convertible_to<Scalar>;
};

} // namespace ctrlpp

#endif
