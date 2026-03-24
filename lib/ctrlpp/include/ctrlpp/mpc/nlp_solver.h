#ifndef HPP_GUARD_CTRLPP_MPC_NLP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_NLP_SOLVER_H

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Core>

#include <concepts>
#include <functional>
#include <span>

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

template <typename S>
concept nlp_solver = requires { typename S::scalar_type; } && requires(S solver, const nlp_problem<typename S::scalar_type>& prob, const nlp_update<typename S::scalar_type>& upd) {
    { solver.setup(prob) } -> std::same_as<void>;
    { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
};

} // namespace ctrlpp

#endif
