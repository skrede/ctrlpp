#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_QP_SOLVER_H

/// @brief Linear-MPC QP backend over argmin's sparse operator-splitting solver.
///
/// This is the argmin-native alternative to `osqp_solver`: it models the same
/// `qp_solver` concept, so `mpc<Scalar, NX, NU, argmin_qp_solver>` compiles and
/// runs without the vendored OSQP C library. Both solve the canonical OSQP-form
/// QP `min ½·xᵀPx + qᵀx  s.t.  l ≤ Ax ≤ u`; argmin ships a header-only C++
/// implementation of the same operator-splitting algorithm (argmin Phases
/// 66/66.1, SEED-042), returning results on a typed error channel rather than
/// through C-pointer ownership.
///
/// The adapter is only defined when argmin's QP header is present. ctrlpp's
/// production argmin pin predates that header, so at the pinned SHA this file is
/// an empty translation unit and no `argmin_qp_solver` symbol exists; a local
/// argmin checkout (or a bumped pin) that ships `argmin/qp/sparse_admm_qp.h`
/// activates it. Promoting it to always-on is a pin-bump decision tracked in the
/// roadmap, not a code change here.
///
/// @cite stellato2020 -- Stellato et al., "OSQP: An Operator Splitting Solver for Quadratic Programs", Math. Prog. Comp. 12(4), 2020

#include "ctrlpp/config.h"

#if defined(CTRLPP_HAS_ARGMIN) && __has_include(<argmin/qp/sparse_admm_qp.h>)

#define CTRLPP_HAS_ARGMIN_QP 1

#include "ctrlpp/expected.h"
#include "ctrlpp/mpc/qp_types.h"

#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/qp/qp_types.h>
#include <argmin/options/sparse_qp_options.h>

#include <cstdint>

namespace ctrlpp
{

/// @brief Structured failure modes for `argmin_qp_solver::try_setup`.
///
///  * pose_failed : argmin rejected the problem at pose time (dimension
///                  mismatch, non-finite data, invalid bounds, or an
///                  unposeable KKT system). No factorization exists, so no
///                  solve can run.
enum class argmin_qp_setup_error : std::uint8_t
{
    pose_failed,
};

/// @brief `qp_solver`-conforming policy wrapping `argmin::sparse_admm_qp_solver`.
///
/// Setup poses and factorizes the problem once (`solve_into`); each subsequent
/// `solve(update)` reuses the frozen Ruiz factors and factorization through
/// argmin's `resolve_into`, warm-starting from the retained iterate exactly as
/// the OSQP backend does. The constructor mirrors `osqp_solver`'s knob roster so
/// the two backends are interchangeable at a call site.
class argmin_qp_solver
{
public:
    using scalar_type = double;

    explicit argmin_qp_solver(double eps_abs = 1e-3, double eps_rel = 1e-3, int max_iter = 4000, bool /*verbose*/ = false, bool warm_starting = true, bool polishing = true)
    {
        // This binds only argmin's stable QP contract: tolerances, iteration
        // budget, and warm-start. Per argmin coordination (argmin-ctrlpp_126-127,
        // SEED-084) the operator-splitting knobs (rho / sigma / alpha /
        // adaptive_rho) are the volatile surface argmin intends to demote to an
        // opt-in sub-struct, so this policy deliberately never touches them --
        // a future reshape of those knobs is a no-op here.
        opts_.eps_abs = eps_abs;
        opts_.eps_rel = eps_rel;
        opts_.max_iterations = static_cast<std::uint16_t>(max_iter);
        opts_.warm_start = warm_starting;
        opts_.polish = polishing;
    }

    /// @brief Preset constructor: `accuracy` polishes, `speed` does not. See
    /// `qp_preset`. For warm-resolve MPC `speed` is the better choice -- the
    /// unpolished iterate already meets the control tolerance at a fraction of
    /// the per-solve cost.
    explicit argmin_qp_solver(qp_preset preset)
        : argmin_qp_solver(1e-3, 1e-3, 4000, false, true, preset == qp_preset::accuracy)
    {
    }

    /// @brief Fallible setup: poses and factorizes the problem once.
    [[nodiscard]] auto try_setup(const qp_problem<double>& problem) -> ctrlpp::expected<void, argmin_qp_setup_error>
    {
        if(auto err = solver_.solve_into(problem.P, problem.q, problem.A, problem.l, problem.u, last_, opts_))
            return ctrlpp::unexpected(argmin_qp_setup_error::pose_failed);
        posed_ = true;
        return {};
    }

    /// @brief Vectors-only resolve reusing the frozen factorization.
    [[nodiscard]] auto solve(const qp_update<double>& update) -> qp_result<double>
    {
        if(!posed_)
            return qp_result<double>{.status = solve_status::error};

        if(update.warm_x.size() > 0 && update.warm_y.size() > 0)
            solver_.warm_start(update.warm_x, update.warm_y);

        if(auto err = solver_.resolve_into(update.q, update.l, update.u, last_, opts_))
            return qp_result<double>{.status = solve_status::error};

        return translate(last_);
    }

private:
    static auto translate(const argmin::qp_result<double>& r) -> qp_result<double>
    {
        qp_result<double> out;
        out.status = translate_status(r.status);
        out.x = r.x;
        out.y = r.y;
        out.objective = r.objective_value;
        out.solve_time = double{0}; // argmin takes no timing measurements
        out.iterations = r.iterations;
        out.primal_residual = r.primal_residual;
        out.dual_residual = r.dual_residual;
        return out;
    }

    static auto translate_status(argmin::qp_solve_status s) -> solve_status
    {
        switch(s)
        {
        case argmin::qp_solve_status::solved:
            return solve_status::optimal;
        case argmin::qp_solve_status::solved_inaccurate:
            return solve_status::solved_inaccurate;
        case argmin::qp_solve_status::max_iterations:
            return solve_status::max_iterations;
        case argmin::qp_solve_status::primal_infeasible:
            return solve_status::infeasible;
        case argmin::qp_solve_status::dual_infeasible:
            return solve_status::unbounded;
        }
        return solve_status::error;
    }

    argmin::sparse_admm_qp_solver<double> solver_{};
    argmin::sparse_qp_options opts_{};
    argmin::qp_result<double> last_{};
    bool posed_{false};
};

}

#endif // CTRLPP_HAS_ARGMIN && __has_include(<argmin/qp/sparse_admm_qp.h>)

#endif
