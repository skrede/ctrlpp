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

template <typename Config>
auto admit_boundary_velocities(Config const& cfg) -> ctrlpp::expected<Config, trajectory_error>
{
    using Scalar = decltype(Config::v_max);
    // k roundings carry a relative error of at most gamma_k = k u / (1 - k u), defined for k u < 1 (@cite higham2002, Lemma 3.1).
    // The admitted k eps = 2 k u is at least gamma_k exactly while k eps <= 1, so a larger count describes no chain of roundings.
    if (!(static_cast<Scalar>(cfg.boundary_velocity_rounding_ops) * std::numeric_limits<Scalar>::epsilon() <= Scalar{1})) {
        return ctrlpp::unexpected(trajectory_error::boundary_velocity_rounding_ops_out_of_range);
    }
    auto const v0 = admit_boundary_velocity(cfg.v0, cfg.v_max, cfg.boundary_velocity_rounding_ops);
    auto const v1 = admit_boundary_velocity(cfg.v1, cfg.v_max, cfg.boundary_velocity_rounding_ops);
    if (!v0.has_value() || !v1.has_value()) {
        return ctrlpp::unexpected(trajectory_error::boundary_velocity_exceeds_limit);
    }
    auto admitted = cfg;
    admitted.v0 = *v0;
    admitted.v1 = *v1;
    return admitted;
}

}

#endif
