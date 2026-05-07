#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_BENCH_SINGLE_SHOOTING_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_BENCH_SINGLE_SHOOTING_H

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Core>

#include <span>
#include <cstddef>
#include <utility>
#include <functional>

namespace bench
{

/// Build a single-shooting NLP for box-constrained solvers.
/// Decision variable: u_0, ..., u_{N-1} (N*NU total).
/// Objective: forward-simulate from x0, accumulate quadratic stage cost plus
/// terminal cost. No equality constraints (continuity is enforced implicitly
/// by the forward simulation embedded in the cost). Bounds are inherited from
/// u_min / u_max applied to every input slot.
template <std::size_t NX, std::size_t NU, typename Dynamics>
inline auto build_single_shooting_problem(
    Dynamics dynamics,
    const Eigen::Matrix<double, NX, 1>& x0,
    int horizon,
    const Eigen::Matrix<double, NX, NX>& Q,
    const Eigen::Matrix<double, NU, NU>& R,
    double u_min,
    double u_max) -> ctrlpp::nlp_problem<double>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    int n_vars = horizon * nu;

    // Cost: forward simulate and accumulate quadratic cost
    auto cost_fn = [=](std::span<const double> z) -> double
    {
        double total = 0.0;
        Eigen::Matrix<double, NX, 1> xk = x0;

        for(int k = 0; k < horizon; ++k)
        {
            Eigen::Map<const Eigen::Matrix<double, NU, 1>> uk(z.data() + k * nu);
            total += 0.5 * xk.dot(Q * xk) + 0.5 * uk.dot(R * uk);
            xk = dynamics(xk, uk);
        }

        // Terminal cost
        total += 0.5 * xk.dot(Q * xk);
        return total;
    };

    std::function<double(std::span<const double>)> cost = cost_fn;

    // Gradient via finite differences
    std::function<void(std::span<const double>, std::span<double>)> gradient =
        [cost](std::span<const double> z, std::span<double> grad)
    {
        ctrlpp::detail::finite_diff_gradient<double>(cost, z, grad);
    };

    // No constraints
    std::function<void(std::span<const double>, std::span<double>)> constraints =
        [](std::span<const double>, std::span<double>) {};

    // Variable bounds: input bounds on all decision variables
    Eigen::VectorXd x_lower = Eigen::VectorXd::Constant(n_vars, u_min);
    Eigen::VectorXd x_upper = Eigen::VectorXd::Constant(n_vars, u_max);

    // No constraint bounds (n_constraints = 0)
    Eigen::VectorXd c_lower;
    Eigen::VectorXd c_upper;

    return ctrlpp::nlp_problem<double>{
        .n_vars = n_vars,
        .n_constraints = 0,
        .cost = std::move(cost),
        .gradient = std::move(gradient),
        .constraints = std::move(constraints),
        .x_lower = std::move(x_lower),
        .x_upper = std::move(x_upper),
        .c_lower = std::move(c_lower),
        .c_upper = std::move(c_upper)};
}

}

#endif
