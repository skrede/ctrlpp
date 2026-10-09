#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_BOUNDARY_VELOCITY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_BOUNDARY_VELOCITY_H

#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <cstdint>

namespace ctrlpp::detail
{

// Only the caller can count the roundings its own chain put on a boundary velocity, so the admitted
// excess is that count in units of eps * v_max rather than a constant chosen here.
template <typename Scalar>
auto admit_boundary_velocity(Scalar v, Scalar v_max, std::uint32_t rounding_ops)
    -> ctrlpp::expected<Scalar, trajectory_error>
{
    if (std::abs(v) <= v_max) {
        return v;
    }
    auto const allowance = static_cast<Scalar>(rounding_ops) * std::numeric_limits<Scalar>::epsilon() * v_max;
    if (std::abs(v) - v_max <= allowance) {
        return std::copysign(v_max, v);
    }
    return ctrlpp::unexpected(trajectory_error::boundary_velocity_exceeds_limit);
}

}

#endif
