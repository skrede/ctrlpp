#ifndef HPP_GUARD_CTRLPP_MPC_NLP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_NLP_SOLVER_H

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_types.h"

#include <Eigen/Core>

#include <span>
#include <cmath>
#include <cstdint>
#include <utility>
#include <variant>
#include <concepts>
#include <functional>

namespace ctrlpp
{

template <typename Scalar>
struct nlp_problem
{
    Eigen::Index n_vars;
    Eigen::Index n_constraints;
    std::function<Scalar(std::span<const Scalar>)> cost_cb;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient_cb;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints_cb;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraint_jacobian_cb;
    Eigen::VectorX<Scalar> x_lower;
    Eigen::VectorX<Scalar> x_upper;
    Eigen::VectorX<Scalar> c_lower;
    Eigen::VectorX<Scalar> c_upper;
};

// Only the decision dimension is compile-time: argmin's fixed-dimension policy already preallocates its
// constraint-axis buffers, so pinning NV alone removes the decision-vector allocations from the solve path.
template <typename Scalar, int NV>
struct nlp_problem_static
{
    static constexpr int problem_dimension = NV;

    Eigen::Index n_vars;
    Eigen::Index n_constraints;
    std::function<Scalar(std::span<const Scalar>)> cost_cb;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient_cb;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints_cb;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraint_jacobian_cb;
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
    nlp_result()
        : status{}
        , x{}
        , objective{}
        , solve_time{}
        , iterations{}
        , primal_residual{}
        , stop_criterion{nlp_stop_criterion::unreported}
    {
    }

    solve_status status;
    Eigen::VectorX<Scalar> x;
    Scalar objective;
    Scalar solve_time;
    std::int32_t iterations;
    Scalar primal_residual;
    nlp_stop_criterion stop_criterion;
};

template <typename S>
concept nlp_solver = requires { typename S::scalar_type; } && requires(S solver, const nlp_update<typename S::scalar_type>& upd) {
    { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
} && requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
    { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
};

template <typename SetupError>
using nmpc_construction_error = std::variant<controller_construction_error, nlp_formulation_error, SetupError>;

namespace detail
{

template <typename Solver, typename Problem>
auto setup_nlp_solver(Solver& solver, const Problem& problem) -> decltype(solver.setup(problem))
{
    return solver.setup(problem);
}

template <typename Solver, typename Problem>
using nlp_setup_error_t = typename decltype(setup_nlp_solver(std::declval<Solver&>(), std::declval<const Problem&>()))::error_type;

// The primal is a separate argument because a write-into solve leaves `result.x` empty and the primal in the
// caller's buffer.
template <typename Scalar>
bool nlp_result_is_usable(const nlp_result<Scalar>& result, const Eigen::VectorX<Scalar>& primal, Eigen::Index primal_size)
{
    return primal.size() >= primal_size && primal.head(primal_size).allFinite() && std::isfinite(result.objective) && std::isfinite(result.solve_time) &&
            std::isfinite(result.primal_residual) && result.iterations >= 0;
}

}

template <typename S>
concept nlp_stepper = nlp_solver<S> && requires(S solver, const nlp_update<typename S::scalar_type>& upd, std::int32_t max_steps) {
    { solver.step(upd, max_steps) } -> std::same_as<nlp_result<typename S::scalar_type>>;
};

}

#endif
