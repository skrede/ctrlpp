#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_STATIC_SOLVE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_STATIC_SOLVE_H

#include "ctrlpp/mpc/detail/nmpc_static_state.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/diagnostics.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, typename Solver, typename Dynamics, std::size_t NC, std::size_t NTC>
class nmpc_static_solve : protected nmpc_static_state<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

    using base = nmpc_static_state<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>;

    // The static (strict-zero) argmin bridge has an allocation-free solve_into; the runtime-erased solvers have only
    // solve(), whose by-value primal allocates.
    static constexpr bool solver_has_solve_into = requires(Solver& s, const nlp_update<Scalar>& u, Eigen::VectorX<Scalar>& out) { s.solve_into(u, out); };

public:
    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0)
    {
        for(auto& ref : this->m_state->x_ref)
            ref.setZero();
        return solve_impl(x0);
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NX>& x_ref)
    {
        for(auto& ref : this->m_state->x_ref)
            ref = x_ref;
        return solve_impl(x0);
    }

    void set_applied_input(const Vector<Scalar, NU>& u) { this->m_u_prev = u; }

    mpc_diagnostics<Scalar> diagnostics() const { return this->m_last_diagnostics; }
    const typename base::problem_type& problem() const { return *this->m_problem; }
    const Eigen::VectorX<Scalar>& last_solution() const { return this->m_last_solution; }

protected:
    using base::horizon;
    using base::problem_dimension;
    using typename base::problem_type;

    using base::base;

private:
    expected<solve_output<Scalar, NU>, solver_error> solve_impl(const Vector<Scalar, NX>& x0)
    {
        if(this->m_setup_failed)
            return unexpected(solver_error::setup_incomplete);

        this->m_state->x0 = x0;
        this->m_state->u_prev = this->m_u_prev;

        this->m_update.x0 = this->m_warm_z;

        const nlp_result<Scalar> result = dispatch_solve();
        populate_diagnostics(result);
        return accept_result(result);
    }

    expected<solve_output<Scalar, NU>, solver_error> accept_result(const nlp_result<Scalar>& result)
    {
        switch(result.status)
        {
            case solve_status::optimal:
                return finish_solve(result, solve_result_status::converged);
            case solve_status::solved_inaccurate:
                return finish_solve(result, solve_result_status::solved_inaccurate);
            case solve_status::max_iterations:
            case solve_status::time_limit:
                return finish_solve(result, solve_result_status::budget_exhausted);
            case solve_status::infeasible:
                return unexpected(solver_error::infeasible);
            case solve_status::unbounded:
            case solve_status::non_convex:
            case solve_status::error:
            default:
                return unexpected(solver_error::invalid_problem);
        }
    }

    // The solve leaves the accepted primal untouched: solve_into writes the pre-sized candidate buffer, and solve()
    // leaves its by-value primal in the result, so a primal of the wrong length resizes no member.
    nlp_result<Scalar> dispatch_solve()
    {
        if constexpr(solver_has_solve_into)
            return this->m_solver.solve_into(this->m_update, this->m_candidate);
        else
            return this->m_solver.solve(this->m_update);
    }

    const Eigen::VectorX<Scalar>& dispatched_primal(const nlp_result<Scalar>& result) const
    {
        if constexpr(solver_has_solve_into)
            return this->m_candidate;
        else
            return result.x;
    }

    // Only an accepted primal is committed: the candidate by a swap, which allocates nothing and keeps the candidate
    // sized for the next write, or the posed prefix of a by-value primal. A refused result touches no member.
    expected<solve_output<Scalar, NU>, solver_error> finish_solve(const nlp_result<Scalar>& result, solve_result_status status)
    {
        if(!detail::nlp_result_is_usable(result, dispatched_primal(result), static_cast<Eigen::Index>(problem_dimension)))
        {
            // The backend's status already reached the diagnostics; correct it, since its answer is discarded.
            this->m_last_diagnostics.status = solve_status::invalid_backend_result;
            return unexpected(solver_error::invalid_backend_result);
        }

        if constexpr(solver_has_solve_into)
            this->m_last_solution.swap(this->m_candidate);
        else
            this->m_last_solution = result.x.head(problem_dimension);
        this->m_has_solution = true;
        shift_warm_start(this->m_last_solution);
        return solve_output<Scalar, NU>{.input = extract_first_input(this->m_last_solution), .status = status};
    }

    void populate_diagnostics(const nlp_result<Scalar>& result)
    {
        this->m_last_diagnostics = mpc_diagnostics<Scalar>{.status = result.status,
                .iterations = result.iterations,
                .solve_time = result.solve_time,
                .cost = result.objective,
                .primal_residual = result.primal_residual,
                .dual_residual = Scalar{0},
                .max_constraint_violation = result.primal_residual,
                .max_path_constraint_violation = Scalar{0},
                .max_terminal_constraint_violation = Scalar{0},
                .total_slack = Scalar{0},
                .stop_criterion = result.stop_criterion};
    }

    // The slack-free formulation has no slack block to shift; the tail input is held.
    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        this->m_warm_z = sol;
        for(std::int32_t k = 0; k < horizon; ++k)
            this->m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
        const Eigen::Index u_offset = (horizon + 1) * nx;
        for(std::int32_t k = 0; k < horizon - 1; ++k)
            this->m_warm_z.segment(u_offset + k * nu, nu) = sol.segment(u_offset + (k + 1) * nu, nu);
    }

    Vector<Scalar, NU> extract_first_input(const Eigen::VectorX<Scalar>& sol)
    {
        const Eigen::Index u_offset = (horizon + 1) * nx;
        Vector<Scalar, NU> u0 = sol.segment(u_offset, nu);
        this->m_u_prev = u0;
        return u0;
    }
};

}

#endif
