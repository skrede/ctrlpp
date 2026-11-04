#ifndef HPP_GUARD_CTRLPP_MHE_MHE_NLP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MHE_MHE_NLP_FORMULATION_H

// @cite diehl2009 -- Diehl et al., "Efficient Numerical Methods for Nonlinear MPC and Moving Horizon Estimation", 2009
// @cite rawlings2017 -- Rawlings, Mayne & Diehl, "Model Predictive Control: Theory, Computation, and Design", 2nd ed., 2017, Ch. 4 (State Estimation)
// @cite rao2003 -- Rao et al., "Constrained State Estimation for Nonlinear Discrete-Time Systems", 2003
// @cite kuhl2011 -- Kuhl et al., 2011 (real-time arrival-cost reassembly)

#include "ctrlpp/mhe/mhe_config.h"

#include "ctrlpp/mhe/detail/nmhe_terms.h"
#include "ctrlpp/mhe/detail/nmhe_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_solver.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"
#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"

#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Dense>

#include <span>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <functional>

namespace ctrlpp
{

// The estimator writes this before each solve; the formulation's closures read it through a shared pointer.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N>
struct nmhe_formulation_state
{
    nmhe_formulation_state()
        : arrival_state{Vector<Scalar, NX>::Zero()}
        , arrival_P_inv{Matrix<Scalar, NX, NX>::Identity()}
        , u_window{}
        , z_window{}
    {
    }

    Vector<Scalar, NX> arrival_state;
    Matrix<Scalar, NX, NX> arrival_P_inv;
    std::array<Vector<Scalar, NU>, N> u_window;
    std::array<Vector<Scalar, NY>, N + 1> z_window;
};

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC>
void set_nmhe_variable_bounds(nlp_problem<Scalar>& problem, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config, const nmhe_layout<NX, N>& layout)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int Ni = static_cast<int>(N);
    problem.x_lower = Eigen::VectorX<Scalar>::Constant(layout.n_vars, -std::numeric_limits<Scalar>::infinity());
    problem.x_upper = Eigen::VectorX<Scalar>::Constant(layout.n_vars, std::numeric_limits<Scalar>::infinity());

    if(config.x_min)
    {
        for(std::int32_t k = 0; k <= Ni; ++k)
            problem.x_lower.segment(layout.x_offset + k * nx, nx) = *config.x_min;
    }
    if(config.x_max)
    {
        for(std::int32_t k = 0; k <= Ni; ++k)
            problem.x_upper.segment(layout.x_offset + k * nx, nx) = *config.x_max;
    }

    if(layout.has_path_slack)
        problem.x_lower.segment(layout.path_slack_offset, layout.n_path_slack).setZero();
}

template <typename Scalar, std::size_t NX, std::size_t N>
void set_nmhe_constraint_bounds(nlp_problem<Scalar>& problem, const nmhe_layout<NX, N>& layout)
{
    problem.c_lower = Eigen::VectorX<Scalar>::Zero(layout.n_constraints);
    problem.c_upper = Eigen::VectorX<Scalar>::Zero(layout.n_constraints);

    for(Eigen::Index i = layout.path_con_start; i < layout.path_con_start + layout.n_path_con; ++i)
        problem.c_lower[i] = -std::numeric_limits<Scalar>::infinity();

    for(Eigen::Index i = layout.residual_con_start; i < layout.residual_con_start + layout.n_residual_con; ++i)
        problem.c_lower[i] = -std::numeric_limits<Scalar>::infinity();
}

// The dynamics enter through the weighted process-noise residual w_k = x_{k+1} - f(x_k, u_k) in the cost, not a
// continuity equality, which would pin every w_k to zero and make the process weighting irrelevant.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC = 0, typename Dynamics, typename Measurement>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
nlp_problem<Scalar> build_nmhe_problem(const Dynamics& dynamics, const Measurement& measurement, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config,
        std::shared_ptr<nmhe_formulation_state<Scalar, NX, NU, NY, N>> state, const Matrix<Scalar, NX, NX>& Q_inv, const Matrix<Scalar, NY, NY>& R_inv)
{
    const nmhe_layout<NX, N> layout = make_nmhe_layout(config);
    std::function<Scalar(std::span<const Scalar>)> cost = nmhe_cost(dynamics, measurement, config, state, Q_inv, R_inv, layout);
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient = [cost, scratch = std::vector<Scalar>(static_cast<std::size_t>(layout.n_vars))](
                                                                                       std::span<const Scalar> z, std::span<Scalar> grad) mutable
    { finite_diff_gradient<Scalar>(cost, z, grad, std::span<Scalar>{scratch.data(), scratch.size()}); };
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints = nmhe_constraints(measurement, config, state, layout);

    // No analytic constraint Jacobian is supplied, so the backend finite-differences it.
    nlp_problem<Scalar> problem{.n_vars = layout.n_vars,
            .n_constraints = layout.n_constraints,
            .cost_cb = std::move(cost),
            .gradient_cb = std::move(gradient),
            .constraints_cb = std::move(constraints),
            .constraint_jacobian_cb = nullptr,
            .x_lower = {},
            .x_upper = {},
            .c_lower = {},
            .c_upper = {}};
    set_nmhe_variable_bounds(problem, config, layout);
    set_nmhe_constraint_bounds(problem, layout);
    return problem;
}

}
}

#endif
