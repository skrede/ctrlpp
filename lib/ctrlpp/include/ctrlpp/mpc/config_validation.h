#ifndef HPP_GUARD_CTRLPP_MPC_CONFIG_VALIDATION_H
#define HPP_GUARD_CTRLPP_MPC_CONFIG_VALIDATION_H

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Core>

#include <cmath>
#include <limits>
#include <optional>

namespace ctrlpp::detail
{

using construction_refusal = std::optional<controller_construction_error>;

template <typename Bound>
auto bounds_are_consistent(const std::optional<Bound>& lower, const std::optional<Bound>& upper) -> bool
{
    using scalar_type = typename Bound::Scalar;
    constexpr scalar_type infinity = std::numeric_limits<scalar_type>::infinity();
    if(lower && (lower->array().isNaN().any() || (lower->array() == infinity).any()))
        return false;
    if(upper && (upper->array().isNaN().any() || (upper->array() == -infinity).any()))
        return false;
    return !(lower && upper) || (lower->array() <= upper->array()).all();
}

template <typename Bound>
auto rate_bound_is_valid(const std::optional<Bound>& du_max) -> bool
{
    return !du_max || (du_max->array() >= 0).all();
}

template <typename Penalty>
auto penalty_is_valid(const Penalty& penalty) -> bool
{
    return penalty.allFinite() && (penalty.array() > 0).all();
}

template <typename Scalar, int N>
auto refused_cost_weight(Eigen::Matrix<Scalar, N, N>& weight) -> construction_refusal
{
    if(!weight.allFinite())
        return controller_construction_error::non_finite_cost_weight;
    // Rounding drift is symmetrized so the posed cost and its analytic gradient
    // share one weight; asymmetry beyond the drift bound is structure and refused.
    if(!symmetrize_within_drift(weight))
        return controller_construction_error::asymmetric_cost_weight;
    if(!psd_decomposition(weight).psd)
        return controller_construction_error::indefinite_cost_weight;
    return std::nullopt;
}

template <typename Config>
auto refused_bounds(const Config& config) -> construction_refusal
{
    if(!bounds_are_consistent(config.u_min, config.u_max))
        return controller_construction_error::invalid_input_bounds;
    if(!bounds_are_consistent(config.x_min, config.x_max))
        return controller_construction_error::invalid_state_bounds;
    if(!rate_bound_is_valid(config.du_max))
        return controller_construction_error::invalid_rate_bound;
    return std::nullopt;
}

template <typename Config>
auto refused_cost_weights(Config& config) -> construction_refusal
{
    if(const auto refusal = refused_cost_weight(config.Q))
        return refusal;
    if(const auto refusal = refused_cost_weight(config.R))
        return refusal;
    if(config.Qf)
        return refused_cost_weight(*config.Qf);
    return std::nullopt;
}

template <typename Config>
auto nmpc_penalties_are_valid(const Config& config) -> bool
{
    const bool path_slack = config.soft_constraints && config.path_penalty.size() > 0 && config.path_constraint;
    const bool terminal_slack = config.soft_constraints && config.terminal_penalty.size() > 0 && config.terminal_constraint;
    return (!path_slack || penalty_is_valid(config.path_penalty)) && (!terminal_slack || penalty_is_valid(config.terminal_penalty));
}

// The uniform penalty is checked even while a per-state one replaces it: it is
// a configuration field with a valid default, so an invalid value is a mistake.
template <typename Config>
auto mpc_penalties_are_valid(const Config& config) -> bool
{
    const bool soft = (config.x_min || config.x_max) && !config.hard_state_constraints;
    const bool uniform = std::isfinite(config.soft_penalty) && config.soft_penalty > 0;
    return !soft || (uniform && (!config.soft_state_penalty || penalty_is_valid(*config.soft_state_penalty)));
}

/// @brief Admit a nonlinear MPC configuration, or name the first rule it
/// violates. Rules, checked in order: input bounds, state bounds, rate bound,
/// then Q, R and a set Qf (each finite, symmetric within rounding drift,
/// positive semidefinite), then the penalties of the active soft constraints.
/// The admitted copy carries the drift-symmetrized weights.
template <typename Config>
auto admit_nmpc_config(const Config& config) -> expected<Config, controller_construction_error>
{
    Config admitted{config};
    if(const auto refusal = refused_bounds(admitted))
        return unexpected(*refusal);
    if(const auto refusal = refused_cost_weights(admitted))
        return unexpected(*refusal);
    if(!nmpc_penalties_are_valid(admitted))
        return unexpected(controller_construction_error::invalid_constraint_penalty);
    return admitted;
}

/// @brief Admit a linear MPC configuration under the same rules and order as
/// `admit_nmpc_config`; the penalties are the slack weights of the soft state
/// constraints.
template <typename Config>
auto admit_mpc_config(const Config& config) -> expected<Config, controller_construction_error>
{
    Config admitted{config};
    if(const auto refusal = refused_bounds(admitted))
        return unexpected(*refusal);
    if(const auto refusal = refused_cost_weights(admitted))
        return unexpected(*refusal);
    if(!mpc_penalties_are_valid(admitted))
        return unexpected(controller_construction_error::invalid_constraint_penalty);
    return admitted;
}

}

#endif
