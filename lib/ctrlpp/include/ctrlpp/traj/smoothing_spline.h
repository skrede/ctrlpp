#ifndef HPP_GUARD_CTRLPP_TRAJ_SMOOTHING_SPLINE_H
#define HPP_GUARD_CTRLPP_TRAJ_SMOOTHING_SPLINE_H

/// @brief Smoothing spline with configurable mu tradeoff parameter.
///
/// Stub -- to be implemented in GREEN phase.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.4.5

#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <vector>

namespace ctrlpp
{

template <typename Scalar>
class smoothing_spline
{
  public:
    struct config
    {
        std::vector<Scalar> times;
        std::vector<Scalar> positions;
        Scalar mu{Scalar{0.5}};
    };

    explicit smoothing_spline(config const& /*cfg*/)
    {
        // Stub -- intentionally not implemented
    }

    auto evaluate(Scalar /*t*/) const -> trajectory_point<Scalar, 1>
    {
        return {};
    }

    auto duration() const -> Scalar { return Scalar{0}; }

  private:
    std::vector<Scalar> times_;
};

static_assert(trajectory_segment<smoothing_spline<double>, double, 1>);

} // namespace ctrlpp

#endif
