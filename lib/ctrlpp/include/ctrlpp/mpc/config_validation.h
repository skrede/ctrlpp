#ifndef HPP_GUARD_CTRLPP_MPC_CONFIG_VALIDATION_H
#define HPP_GUARD_CTRLPP_MPC_CONFIG_VALIDATION_H

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Core>

#include <cmath>
#include <limits>
#include <cstddef>
#include <variant>
#include <optional>

namespace ctrlpp::detail
{

using construction_refusal = std::optional<controller_construction_error>;

template <typename Bound>
bool bounds_are_consistent(const std::optional<Bound>& lower, const std::optional<Bound>& upper)
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
bool rate_bound_is_valid(const std::optional<Bound>& du_max)
{
    return !du_max || (du_max->array() >= 0).all();
}

template <typename Penalty>
bool penalty_is_valid(const Penalty& penalty)
{
    return penalty.allFinite() && (penalty.array() > 0).all();
}

struct weight_refusals
{
    controller_construction_error non_finite;
    controller_construction_error asymmetric;
    controller_construction_error indefinite;
};

inline constexpr weight_refusals cost_weight_refusals{
        controller_construction_error::non_finite_cost_weight, controller_construction_error::asymmetric_cost_weight, controller_construction_error::indefinite_cost_weight};

inline constexpr weight_refusals terminal_set_refusals{
        controller_construction_error::non_finite_terminal_set, controller_construction_error::asymmetric_terminal_set, controller_construction_error::indefinite_terminal_set};

template <typename Scalar, int N>
construction_refusal refused_weight(Eigen::Matrix<Scalar, N, N>& weight, const weight_refusals& refusals)
{
    if(!weight.allFinite())
        return refusals.non_finite;
    // Rounding drift is symmetrized so every triangle the posed problem reads is
    // the same matrix; asymmetry beyond the drift bound is structure and refused.
    if(!symmetrize_within_drift(weight))
        return refusals.asymmetric;
    if(!psd_decomposition(weight).psd)
        return refusals.indefinite;
    return std::nullopt;
}

template <typename Config>
construction_refusal refused_bounds(const Config& config)
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
construction_refusal refused_cost_weights(Config& config)
{
    if(const auto refusal = refused_weight(config.Q, cost_weight_refusals))
        return refusal;
    if(const auto refusal = refused_weight(config.R, cost_weight_refusals))
        return refusal;
    if(config.Qf)
        return refused_weight(*config.Qf, cost_weight_refusals);
    return std::nullopt;
}

template <typename Scalar, std::size_t NX>
construction_refusal refused_terminal_set(ellipsoidal_set<Scalar, NX>& set)
{
    if(const auto refusal = refused_weight(set.P, terminal_set_refusals))
        return refusal;
    if(!std::isfinite(set.alpha))
        return controller_construction_error::non_finite_terminal_set;
    if(set.alpha < Scalar{0})
        return controller_construction_error::negative_terminal_set_level;
    return std::nullopt;
}

template <typename Scalar, std::size_t NX>
construction_refusal refused_terminal_set(const polytopic_set<Scalar, NX>& set)
{
    if(set.H.rows() != set.h.size())
        return controller_construction_error::mismatched_terminal_set;
    if(!set.H.allFinite() || set.h.array().isNaN().any())
        return controller_construction_error::non_finite_terminal_set;
    return std::nullopt;
}

template <typename Config>
bool nmpc_penalties_are_valid(const Config& config)
{
    const bool path_slack = config.soft_constraints && config.path_penalty.size() > 0 && config.path_constraint_cb;
    const bool terminal_slack = config.soft_constraints && config.terminal_penalty.size() > 0 && config.terminal_constraint_cb;
    return (!path_slack || penalty_is_valid(config.path_penalty)) && (!terminal_slack || penalty_is_valid(config.terminal_penalty));
}

// The uniform penalty is checked even while a per-state one replaces it: it is
// a configuration field with a valid default, so an invalid value is a mistake.
template <typename Config>
bool mpc_penalties_are_valid(const Config& config)
{
    const bool soft = (config.x_min || config.x_max) && !config.hard_state_constraints;
    const bool uniform = std::isfinite(config.soft_penalty) && config.soft_penalty > 0;
    return !soft || (uniform && (!config.soft_state_penalty || penalty_is_valid(*config.soft_state_penalty)));
}

// The first violated rule is returned, in the order checked here; the admitted copy carries the drift-symmetrized
// weights.
template <typename Config>
expected<Config, controller_construction_error> admit_nmpc_config(const Config& config)
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

// The nonlinear admission's rules and order, then the terminal set; the admitted copy carries the symmetrized P.
template <typename Config>
expected<Config, controller_construction_error> admit_mpc_config(const Config& config)
{
    Config admitted{config};
    if(const auto refusal = refused_bounds(admitted))
        return unexpected(*refusal);
    if(const auto refusal = refused_cost_weights(admitted))
        return unexpected(*refusal);
    if(!mpc_penalties_are_valid(admitted))
        return unexpected(controller_construction_error::invalid_constraint_penalty);
    if(admitted.terminal_constraint_set)
        if(const auto refusal = std::visit([](auto& set) -> construction_refusal { return refused_terminal_set(set); }, *admitted.terminal_constraint_set))
            return unexpected(*refusal);
    return admitted;
}

}

#endif
