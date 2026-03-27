#ifndef HPP_GUARD_CTRLPP_TRAJ_PIECEWISE_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJ_PIECEWISE_TRAJECTORY_H

/// @brief Variadic piecewise trajectory composition.
///
/// Composes heterogeneous trajectory segments into a single trajectory that
/// evaluates seamlessly across segment boundaries with automatic time offset.
/// The piecewise_trajectory itself satisfies trajectory_segment, enabling
/// recursive composition.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Ch. 4 (composite trajectories)

#include "ctrlpp/traj/cubic_trajectory.h"
#include "ctrlpp/traj/detail/piecewise_dispatch.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <cstddef>

namespace ctrlpp
{

/// @brief Composes heterogeneous trajectory segments into a single trajectory.
///
/// Segments are evaluated in sequence with automatic time offset. At breakpoints,
/// the next segment is evaluated at local t=0 to ensure continuity.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Sec. 4.1 -- multi-segment trajectory with automatic time offsets
template <typename Scalar, std::size_t ND, typename... Segments>
    requires(trajectory_segment<Segments, Scalar, ND> && ...)
class piecewise_trajectory
    : public detail::piecewise_base<piecewise_trajectory<Scalar, ND, Segments...>, Scalar,
                                    trajectory_point<Scalar, ND>, Segments...>
{
    using base =
        detail::piecewise_base<piecewise_trajectory<Scalar, ND, Segments...>, Scalar,
                               trajectory_point<Scalar, ND>, Segments...>;

public:
    using base::base;
};

/// Verify piecewise_trajectory itself satisfies trajectory_segment (recursive composition).
static_assert(trajectory_segment<
    piecewise_trajectory<double, 1, cubic_trajectory<double, 1>, cubic_trajectory<double, 1>>,
    double, 1>);

} // namespace ctrlpp

#endif
