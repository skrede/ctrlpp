#ifndef HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H

#include "ctrlpp/mpc/qp_types.h"

#include <utility>
#include <concepts>

namespace ctrlpp
{

// @cite stellato2020 -- Stellato et al., "OSQP: An Operator Splitting Solver for Quadratic Programs", Math. Prog. Comp. 12(4), 2020 (default backend that satisfies this concept)
template <typename S>
concept qp_solver = requires { typename S::scalar_type; } && requires(S solver, const qp_update<typename S::scalar_type>& upd) {
    { solver.solve(upd) } -> std::same_as<qp_result<typename S::scalar_type>>;
} && requires(S solver, const qp_problem<typename S::scalar_type>& prob) {
    { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
};

namespace detail
{

template <typename Solver, typename Problem>
auto setup_qp_solver(Solver& solver, const Problem& problem) -> decltype(solver.setup(problem))
{
    return solver.setup(problem);
}

template <typename Solver, typename Problem>
using qp_setup_error_t = typename decltype(setup_qp_solver(std::declval<Solver&>(), std::declval<const Problem&>()))::error_type;

}

}

#endif
