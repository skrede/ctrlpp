#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_PIECEWISE_PATH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_PIECEWISE_PATH_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2009, Ch. 4 (composite trajectories)

#include "ctrlpp/trajectory/path_segment.h"
#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/trajectory/detail/piecewise_dispatch.h"

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 4.1, eq. (4.1) -- multi-segment path construction
template <typename Scalar, typename... Segments>
    requires(path_segment<Segments, Scalar> && ...)
class piecewise_path : public detail::piecewise_base<piecewise_path<Scalar, Segments...>, Scalar, path_point<Scalar>, Segments...>
{
    using base = detail::piecewise_base<piecewise_path<Scalar, Segments...>, Scalar, path_point<Scalar>, Segments...>;

public:
    using base::base;
};

}

#endif
