#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_MPC_SOLVE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_MPC_SOLVE_H

#include "ctrlpp/mpc/detail/mpc_state.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/qp_formulation.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, typename Solver, std::size_t NY, typename Config>
class mpc_solve : protected mpc_state<Scalar, NX, NU, Solver, NY, Config>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

    using base = mpc_state<Scalar, NX, NU, Solver, NY, Config>;

protected:
    using base::base;

    expected<solve_output<Scalar, NU>, solver_error> solve_impl(const Vector<Scalar, NX>& x0)
    {
        rebuild_bounds(x0);
        apply_rate_constraint_update();
        apply_warm_start();
        auto result = this->m_solver.solve(this->m_update);
        populate_diagnostics(result);

        switch(result.status)
        {
            case solve_status::optimal:
                return extract_solution(result, solve_result_status::converged);
            case solve_status::solved_inaccurate:
                return extract_solution(result, solve_result_status::solved_inaccurate);
            case solve_status::max_iterations:
            case solve_status::time_limit:
                return extract_solution(result, solve_result_status::budget_exhausted);
            case solve_status::infeasible:
                return unexpected(solver_error::infeasible);
            case solve_status::unbounded:
            case solve_status::non_convex:
            case solve_status::error:
            default:
                return unexpected(solver_error::invalid_problem);
        }
    }

private:
    void rebuild_bounds(const Vector<Scalar, NX>& x0)
    {
        std::int32_t N = this->m_config.horizon;
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(
                N, x0, this->m_config.x_min, this->m_config.x_max, this->m_config.u_min, this->m_config.u_max, this->m_config.du_max, this->m_config.terminal_constraint_set);
        this->m_update.l = std::move(l);
        this->m_update.u = std::move(u);
    }

    void apply_rate_constraint_update()
    {
        if(!this->m_config.du_max.has_value())
            return;

        // The offset comes from the same derivation that sized the bounds vectors
        // this function writes into, so it cannot address a different block.
        this->m_update.l.segment(this->m_layout.rate_row, nu) = -(*this->m_config.du_max) + this->m_u_prev;
        this->m_update.u.segment(this->m_layout.rate_row, nu) = *this->m_config.du_max + this->m_u_prev;
    }

    void apply_warm_start()
    {
        if(this->m_has_solution)
        {
            this->m_update.warm_x = this->m_last_primal;
            this->m_update.warm_y = this->m_last_dual;
        }
    }

    static mpc_diagnostics<Scalar> diagnostics_of(const qp_result<Scalar>& result)
    {
        mpc_diagnostics<Scalar> d;
        d.status = result.status;
        d.iterations = result.iterations;
        d.solve_time = result.solve_time;
        d.cost = result.objective;
        d.primal_residual = result.primal_residual;
        d.dual_residual = result.dual_residual;
        d.max_constraint_violation = Scalar{0};
        return d;
    }

    void populate_diagnostics(const qp_result<Scalar>& result) { this->m_last_diagnostics = diagnostics_of(result); }

    // A status is not a shape: an accepted status can come with a primal or dual shorter than the posed problem, or with
    // non-finite values, so both are checked before either vector is read, and a refusal leaves the command and the
    // warm start untouched.
    expected<solve_output<Scalar, NU>, solver_error> extract_solution(qp_result<Scalar>& result, solve_result_status status)
    {
        auto const primal_size = this->m_n_dec;
        auto const dual_size = this->m_layout.n_con;
        if(result.x.size() < primal_size || result.y.size() < dual_size || !result.x.head(primal_size).allFinite() || !result.y.head(dual_size).allFinite() ||
                !std::isfinite(result.objective) || !std::isfinite(result.solve_time) || !std::isfinite(result.primal_residual) || !std::isfinite(result.dual_residual) ||
                result.iterations < 0)
        {
            // Correct the diagnostics the backend's own status just populated, so a reader there is not told the
            // solve was optimal when its answer was discarded.
            this->m_last_diagnostics.status = solve_status::invalid_backend_result;
            return unexpected(solver_error::invalid_backend_result);
        }

        this->m_last_primal = std::move(result.x);
        this->m_last_dual = std::move(result.y);
        this->m_has_solution = true;
        measure_constraints();

        Vector<Scalar, NU> u0 = this->m_last_primal.segment(this->m_n_x_total, nu);
        this->m_u_prev = u0;
        return solve_output<Scalar, NU>{.input = std::move(u0), .status = status};
    }

    // The overall figure takes the rows as posed, where a soft state row carries its slack column, so an overshoot the
    // slack absorbs is not counted there; the path figure takes the state bounds as stated, so it is.
    void measure_constraints()
    {
        auto const& lay = this->m_layout;
        auto const& z = this->m_last_primal;
        auto& Ax = this->m_constraint_product;
        Ax.noalias() = this->m_constraint_matrix * z.head(this->m_n_dec);

        auto& d = this->m_last_diagnostics;
        d.max_constraint_violation = bound_violation(Ax, 0);
        d.max_path_constraint_violation = std::max(bound_violation(z.segment(nx, lay.input_row - lay.state_row), lay.state_row),
                bound_violation(Ax.segment(lay.input_row, lay.terminal_row - lay.input_row), lay.input_row));
        d.max_terminal_constraint_violation = bound_violation(Ax.segment(lay.terminal_row, lay.n_con - lay.terminal_row), lay.terminal_row);

        auto const slack_offset = this->m_n_x_total + this->m_config.horizon * nu;
        d.total_slack = z.segment(slack_offset, this->m_n_dec - slack_offset).cwiseAbs().sum();
    }

    template <typename Values>
    Scalar bound_violation(const Eigen::MatrixBase<Values>& values, Eigen::Index row) const
    {
        const Eigen::Index n = values.size();
        if(n == 0)
            return Scalar{0};
        auto const l = this->m_update.l.segment(row, n);
        auto const u = this->m_update.u.segment(row, n);
        return (l - values).cwiseMax(values - u).cwiseMax(Scalar{0}).maxCoeff();
    }
};

}

#endif
