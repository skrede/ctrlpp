#ifndef HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H

// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017, Ch. 8 (NMPC, multiple shooting)
// @cite diehl2002 -- Diehl, Bock, Schloder et al., "Real-Time Optimization and Nonlinear Model Predictive Control of Processes Governed by DAEs", J. Process Control 12(4), 2002 (real-time iteration scheme)

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/detail/nmpc_bounds.h"
#include "ctrlpp/mpc/detail/nmpc_layout.h"
#include "ctrlpp/mpc/detail/nmpc_gradient.h"
#include "ctrlpp/mpc/detail/nmpc_jacobian.h"
#include "ctrlpp/mpc/detail/nmpc_cost_terms.h"
#include "ctrlpp/mpc/detail/nmpc_constraint_terms.h"

#include <Eigen/Dense>

#include <span>
#include <memory>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <functional>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU>
struct nmpc_formulation_state
{
    nmpc_formulation_state()
        : x0{Vector<Scalar, NX>::Zero()}
        , x_ref{}
        , u_prev{Vector<Scalar, NU>::Zero()}
    {
    }

    Vector<Scalar, NX> x0;
    std::vector<Vector<Scalar, NX>> x_ref;
    Vector<Scalar, NU> u_prev;
};

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC = 0, std::size_t NTC = 0, typename Dynamics>
nlp_problem<Scalar> build_nmpc_problem(const Dynamics& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> state)
{
    const nmpc_layout layout = make_nmpc_layout(config);
    std::function<Scalar(std::span<const Scalar>)> cost = nmpc_cost(config, state, layout);
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient = make_nmpc_gradient(config, state, layout, cost);
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints = nmpc_constraints(dynamics, config, state, layout);
    nmpc_bounds<Scalar> x_bounds = make_nmpc_variable_bounds(config, layout);
    nmpc_bounds<Scalar> c_bounds = make_nmpc_constraint_bounds<Scalar>(layout);
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraint_jacobian = nmpc_constraint_jacobian(config, layout, constraints);

    return nlp_problem<Scalar>{.n_vars = layout.n_vars,
            .n_constraints = layout.n_constraints,
            .cost_cb = std::move(cost),
            .gradient_cb = std::move(gradient),
            .constraints_cb = std::move(constraints),
            .constraint_jacobian_cb = std::move(constraint_jacobian),
            .x_lower = std::move(x_bounds.lower),
            .x_upper = std::move(x_bounds.upper),
            .c_lower = std::move(c_bounds.lower),
            .c_upper = std::move(c_bounds.upper)};
}

// Shared by the formulation and the controller's factory, so both refuse the same configurations, horizon first.
template <std::size_t NH, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC>
ctrlpp::expected<void, nlp_formulation_error> static_formulation_refusal(const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
{
    if(config.horizon != static_cast<std::int32_t>(NH))
        return ctrlpp::unexpected(nlp_formulation_error::horizon_mismatch);
    const bool has_path_slack = config.soft_constraints && NC > 0 && config.path_constraint_cb.has_value();
    const bool has_term_slack = config.soft_constraints && NTC > 0 && config.terminal_constraint_cb.has_value();
    if(has_path_slack || has_term_slack)
        return ctrlpp::unexpected(nlp_formulation_error::slack_not_supported);
    return {};
}

// Poses only the slack-free configuration, where n_vars == NV; the constraint count stays runtime.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, std::size_t NC = 0, std::size_t NTC = 0, typename Dynamics>
ctrlpp::expected<nlp_problem_static<Scalar, static_cast<int>((NH + 1) * NX + NH * NU)>, nlp_formulation_error> build_nmpc_problem_static(
        const Dynamics& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> state)
{
    static_assert(NH > 0,
            "Horizon NH must be positive: a zero horizon leaves no input block in the decision vector, so the (NH+1)*NX + NH*NU dimension collapses onto the initial state alone and "
            "every input offset the assembled problem writes lands outside it");

    static constexpr int NV = static_cast<int>((NH + 1) * NX + NH * NU);

    // Refused before build_nmpc_problem runs, since everything it builds is sized from config.horizon.
    if(auto refused = static_formulation_refusal<NH>(config); !refused.has_value())
        return ctrlpp::unexpected(refused.error());

    auto dyn = build_nmpc_problem<Scalar, NX, NU, NC, NTC>(dynamics, config, state);

    nlp_problem_static<Scalar, NV> prob;
    prob.n_vars = dyn.n_vars;
    prob.n_constraints = dyn.n_constraints;
    prob.cost_cb = std::move(dyn.cost_cb);
    prob.gradient_cb = std::move(dyn.gradient_cb);
    prob.constraints_cb = std::move(dyn.constraints_cb);
    prob.constraint_jacobian_cb = std::move(dyn.constraint_jacobian_cb);
    prob.x_lower = dyn.x_lower;
    prob.x_upper = dyn.x_upper;
    prob.c_lower = std::move(dyn.c_lower);
    prob.c_upper = std::move(dyn.c_upper);
    return prob;
}

}
}

#endif
