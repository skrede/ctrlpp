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

/// @brief Compile-time-dimension NLP contract, parallel to the runtime-erased
/// nlp_problem<Scalar>. It carries the decision dimension NV as a static
/// template parameter so the argmin bridge and argmin's compile-time-N solver
/// can size their decision-vector storage with fixed-size Eigen types (the
/// allocation-free static path, Route A / SEED-002).
///
/// Only the DECISION dimension is compile-time here: the bound vectors
/// x_lower / x_upper are fixed-size Eigen::Vector<Scalar, NV>, while the
/// constraint count stays runtime (c_lower / c_upper remain dynamic and
/// n_constraints is an int). argmin's fixed-N NW-SQP policy already keeps its
/// constraint-axis buffers dynamic-but-preallocated, so pinning NV alone
/// removes the decision-vector allocations that dominate the RT hot path. The
/// cost/gradient/constraints/constraint_jacobian callables keep the same
/// std::span shape as nlp_problem (a std::function call does not allocate).
template <typename Scalar, int NV>
struct nlp_problem_static
{
    static constexpr int problem_dimension = NV;

    int n_vars;
    int n_constraints;
    std::function<Scalar(std::span<const Scalar>)> cost;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraint_jacobian;
    Eigen::Vector<Scalar, NV> x_lower;
    Eigen::Vector<Scalar, NV> x_upper;
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
/// NMHE. A solver models the concept with one setup shape: a fallible
/// `setup(problem)` whose result exposes `has_value()` (a
/// `ctrlpp::expected<void, E>` over the backend's own setup-error enum). A
/// backend with no setup failure mode writes a trivially succeeding fallible
/// setup rather than an infallible one, so every backend reports a setup
/// outcome through the same typed channel and none reports it out of band. The
/// per-iteration `solve` always returns an `nlp_result` carrying a plain
/// `solve_status` value.
template <typename S>
concept nlp_solver = requires { typename S::scalar_type; }
    && requires(S solver, const nlp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    }
    && requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
        { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
    };

namespace detail
{

/// @brief Dispatch NLP solver setup. The fallible `setup` is the single shape
/// the concept accepts, so this forwards both the call and the backend's own
/// typed error: the return type is whatever `Solver::setup` returns, a
/// `ctrlpp::expected<void, E>` over that backend's setup-error enum. Nothing is
/// flattened here, so a caller that wants the cause of a setup failure can read
/// it at this seam. Generic over the problem type, so a solver bound to the
/// compile-time-dimension `nlp_problem_static` uses the same helper as one
/// bound to the runtime-erased `nlp_problem`.
template <typename Solver, typename Problem>
auto setup_nlp_solver(Solver& solver, const Problem& problem) -> decltype(solver.setup(problem))
{
    return solver.setup(problem);
}

}

template <typename S>
concept nlp_stepper = nlp_solver<S> &&
    requires(S solver, const nlp_update<typename S::scalar_type>& upd, int max_steps) {
        { solver.step(upd, max_steps) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    };

}

#endif
