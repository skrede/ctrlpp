#ifndef HPP_GUARD_CTRLPP_TRAJ_PIECEWISE_PATH_H
#define HPP_GUARD_CTRLPP_TRAJ_PIECEWISE_PATH_H

/// @brief Variadic piecewise path composition.
///
/// Composes heterogeneous path segments into a single path that evaluates
/// seamlessly across segment boundaries with automatic time offset.
/// The piecewise_path itself satisfies path_segment, enabling recursive composition.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Ch. 4 (composite trajectories)

#include "ctrlpp/traj/detail/piecewise_dispatch.h"
#include "ctrlpp/traj/path_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

namespace ctrlpp
{

/// @brief Composes heterogeneous path segments into a single path.
///
/// Segments are evaluated in sequence with automatic time offset. At breakpoints,
/// the next segment is evaluated at local t=0 to ensure continuity.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Sec. 4.1, eq. (4.1) -- multi-segment path construction
template <typename Scalar, typename... Segments>
    requires(path_segment<Segments, Scalar> && ...)
class piecewise_path
    : public detail::piecewise_base<piecewise_path<Scalar, Segments...>, Scalar, path_point<Scalar>,
                                    Segments...>
{
    using base =
        detail::piecewise_base<piecewise_path<Scalar, Segments...>, Scalar, path_point<Scalar>,
                               Segments...>;

public:
    using base::base;
};

} // namespace ctrlpp

#endif
