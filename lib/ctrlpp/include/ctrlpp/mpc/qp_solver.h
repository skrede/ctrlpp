#ifndef HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H

/// @brief Concept defining the QP solver interface used by linear MPC.
///
/// A solver models the concept with either setup shape: a `void setup(problem)`
/// (infallible or throwing) or a fallible `try_setup(problem)` whose result
/// exposes `has_value()` (an `expected<void, E>`). The per-iteration `solve`
/// always returns a `qp_result` carrying a plain `solve_status` value.
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
    && (requires(S solver, const qp_problem<typename S::scalar_type>& prob) {
            { solver.setup(prob) } -> std::same_as<void>;
        }
        || requires(S solver, const qp_problem<typename S::scalar_type>& prob) {
            { solver.try_setup(prob).has_value() } -> std::convertible_to<bool>;
        });

namespace detail
{

/// @brief Dispatch QP solver setup through `try_setup` when the solver provides
/// it. Returns true when setup succeeded; solvers exposing only the classic
/// `setup(problem)` report success unconditionally (their failures, if any,
/// escape as exceptions).
template <typename Solver, typename Scalar>
[[nodiscard]] auto setup_qp_solver(Solver& solver, const qp_problem<Scalar>& problem) -> bool
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

}

#endif
