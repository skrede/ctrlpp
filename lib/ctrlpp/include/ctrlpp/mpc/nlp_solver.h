#ifndef HPP_GUARD_CTRLPP_MPC_NLP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_NLP_SOLVER_H

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Core>

#include <span>
#include <concepts>
#include <functional>

namespace ctrlpp
{

template <typename Scalar>
struct nlp_problem
{
    int n_vars;
    int n_constraints;
    std::function<Scalar(std::span<const Scalar>)> cost;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraint_jacobian;
    Eigen::VectorX<Scalar> x_lower;
    Eigen::VectorX<Scalar> x_upper;
    Eigen::VectorX<Scalar> c_lower;
    Eigen::VectorX<Scalar> c_upper;
};

template <typename Scalar>
struct nlp_update
{
    Eigen::VectorX<Scalar> x0;
};

template <typename Scalar>
struct nlp_result
{
    solve_status status;
    Eigen::VectorX<Scalar> x;
    Scalar objective;
    Scalar solve_time;
    int iterations;
    Scalar primal_residual;
};

/// @brief Concept defining the NLP solver interface used by nonlinear MPC and
/// NMHE. A solver models the concept with either setup shape: a
/// `void setup(problem)` (infallible or throwing) or a fallible
/// `try_setup(problem)` whose result exposes `has_value()` (an
/// `expected<void, E>`). The per-iteration `solve` always returns an
/// `nlp_result` carrying a plain `solve_status` value.
template <typename S>
concept nlp_solver = requires { typename S::scalar_type; }
    && requires(S solver, const nlp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    }
    && (requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
            { solver.setup(prob) } -> std::same_as<void>;
        }
        || requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
            { solver.try_setup(prob).has_value() } -> std::convertible_to<bool>;
        });

namespace detail
{

/// @brief Dispatch NLP solver setup through `try_setup` when the solver
/// provides it. Returns true when setup succeeded; solvers exposing only the
/// classic `setup(problem)` report success unconditionally (their failures, if
/// any, escape as exceptions).
template <typename Solver, typename Scalar>
[[nodiscard]] auto setup_nlp_solver(Solver& solver, const nlp_problem<Scalar>& problem) -> bool
{
    if constexpr(requires { solver.try_setup(problem); })
    {
        return solver.try_setup(problem).has_value();
    }
    else
    {
        solver.setup(problem);
        return true;
    }
}

}

template <typename S>
concept nlp_stepper = nlp_solver<S> &&
    requires(S solver, const nlp_update<typename S::scalar_type>& upd, int max_steps) {
        { solver.step(upd, max_steps) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    };

}

#endif
