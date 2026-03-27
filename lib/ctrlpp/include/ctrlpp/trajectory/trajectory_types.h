#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_TYPES_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_TYPES_H

/// @brief Core output types for trajectory generation.
///
/// trajectory_point<Scalar, ND> holds position, velocity, and acceleration vectors
/// for an ND-dimensional trajectory. path_point<Scalar> holds normalized scalar
/// values (q, dq, ddq, dddq) for a path mapping tau in [0,1].
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 5.2.1, eq. (2.16)

#include "ctrlpp/types.h"

#include <cstddef>

namespace ctrlpp
{

/// @brief Point on an ND-dimensional trajectory with position, velocity, acceleration.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Sec. 2.1, eq. (2.1)-(2.3) -- trajectory representation
///   with position q(t), velocity q'(t), acceleration q''(t)
template <typename Scalar, std::size_t ND>
struct trajectory_point
{
    Vector<Scalar, ND> position{};
    Vector<Scalar, ND> velocity{};
    Vector<Scalar, ND> acceleration{};
};

/// @brief Normalized path output for tau in [0,1].
///
/// Fields represent normalized position, velocity, acceleration, and jerk.
/// Physical values are obtained via scaling: vel = h/T * dq, acc = h/T^2 * ddq, etc.
template <typename Scalar>
struct path_point
{
    Scalar q{};
    Scalar dq{};
    Scalar ddq{};
    Scalar dddq{};
};

} // namespace ctrlpp

#endif
