#ifndef HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H

/// @brief Concept defining the QP solver interface used by linear MPC.
///
/// A solver models the concept with one setup shape: a fallible
/// `setup(problem)` whose result exposes `has_value()` (a
/// `ctrlpp::expected<void, E>` over the backend's own setup-error enum). A
/// backend with no setup failure mode writes a trivially succeeding fallible
/// setup rather than an infallible one, so every backend reports a setup
/// outcome through the same typed channel and none reports it out of band. The
/// per-iteration `solve` always returns a `qp_result` carrying a plain
/// `solve_status` value.
///
/// @cite stellato2020 -- Stellato et al., "OSQP: An Operator Splitting Solver for Quadratic Programs", Math. Prog. Comp. 12(4), 2020 (default backend that satisfies this concept)

#include "ctrlpp/mpc/qp_types.h"

#include <concepts>

namespace ctrlpp
{

template <typename S>
concept qp_solver = requires { typename S::scalar_type; }
    && requires(S solver, const qp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<qp_result<typename S::scalar_type>>;
    }
    && requires(S solver, const qp_problem<typename S::scalar_type>& prob) {
        { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
    };

namespace detail
{

/// @brief Dispatch QP solver setup. The fallible `setup` is the single shape
/// the concept accepts, so this forwards both the call and the backend's own
/// typed error: the return type is whatever `Solver::setup` returns, a
/// `ctrlpp::expected<void, E>` over that backend's setup-error enum. Nothing is
/// flattened here, so a caller that wants the cause of a setup failure can read
/// it at this seam. Generic over the problem type, so a solver bound to a
/// compile-time-dimension problem uses the same helper as a runtime-erased one.
template <typename Solver, typename Problem>
auto setup_qp_solver(Solver& solver, const Problem& problem) -> decltype(solver.setup(problem))
{
    return solver.setup(problem);
}

}

}

#endif
