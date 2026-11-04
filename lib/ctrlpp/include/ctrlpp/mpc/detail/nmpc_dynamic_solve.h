#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_DYNAMIC_SOLVE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_DYNAMIC_SOLVE_H

#include "ctrlpp/mpc/detail/nmpc_dynamic_result.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/diagnostics.h"

#include <Eigen/Dense>

#include <span>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, typename Solver, typename Dynamics, std::size_t NC, std::size_t NTC>
class nmpc_dynamic_solve : protected nmpc_dynamic_result<Scalar, NX, NU, Solver, Dynamics, NC, NTC>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

    using base = nmpc_dynamic_result<Scalar, NX, NU, Solver, Dynamics, NC, NTC>;

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

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NX>> x_ref)
    {
        if(x_ref.empty())
            return unexpected(solver_error::invalid_problem);
        const auto len = std::min(x_ref.size(), this->m_state->x_ref.size());
        for(std::size_t k = 0; k < len; ++k)
            this->m_state->x_ref[k] = x_ref[k];
        for(std::size_t k = len; k < this->m_state->x_ref.size(); ++k)
            this->m_state->x_ref[k] = x_ref.back();
        return solve_impl(x0);
    }

    void set_applied_input(const Vector<Scalar, NU>& u) { this->m_u_prev = u; }

    expected<std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>, solver_error> trajectory() const
    {
        if(!this->m_has_solution)
            return unexpected(solver_error::setup_incomplete);

        // No length check: m_has_solution is set only after the stored primal was checked to cover the problem
        // dimension these offsets derive from.
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(this->m_N + 1));
        inputs.reserve(static_cast<std::size_t>(this->m_N));

        for(std::int32_t k = 0; k <= this->m_N; ++k)
            states.push_back(this->m_last_solution.segment(k * nx, nx));
        const Eigen::Index u_offset = (this->m_N + 1) * nx;
        for(std::int32_t k = 0; k < this->m_N; ++k)
            inputs.push_back(this->m_last_solution.segment(u_offset + k * nu, nu));

        return std::pair{std::move(states), std::move(inputs)};
    }

    mpc_diagnostics<Scalar> diagnostics() const { return this->m_last_diagnostics; }
    const nlp_problem<Scalar>& problem() const { return *this->m_problem; }
    const Eigen::VectorX<Scalar>& last_solution() const { return this->m_last_solution; }

protected:
    using base::base;

private:
    expected<solve_output<Scalar, NU>, solver_error> solve_impl(const Vector<Scalar, NX>& x0)
    {
        if(this->m_setup_failed)
            return unexpected(solver_error::setup_incomplete);

        this->m_state->x0 = x0;
        this->m_state->u_prev = this->m_u_prev;

        nlp_update<Scalar> update;
        update.x0 = this->m_warm_z;

        auto result = this->m_solver.solve(update);
        this->populate_diagnostics(result);
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

    // A refused result touches no member, so the warm start, the previous input and the trajectory keep the last
    // accepted solve.
    expected<solve_output<Scalar, NU>, solver_error> finish_solve(const nlp_result<Scalar>& result, solve_result_status status)
    {
        const Eigen::VectorX<Scalar>& z = result.x;
        if(!detail::nlp_result_is_usable(result, z, this->m_num_vars))
        {
            // The backend's status already reached the diagnostics; correct it, since its answer is discarded.
            this->m_last_diagnostics.status = solve_status::invalid_backend_result;
            return unexpected(solver_error::invalid_backend_result);
        }

        this->m_last_solution = z;
        this->m_has_solution = true;
        this->populate_constraint_diagnostics(z);
        this->shift_warm_start(z);
        return solve_output<Scalar, NU>{.input = this->extract_first_input(z), .status = status};
    }
};

}

#endif
