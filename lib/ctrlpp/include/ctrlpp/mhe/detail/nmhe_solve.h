#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_SOLVE_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_SOLVE_H

#include "ctrlpp/mhe/detail/nmhe_state.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC, typename Solver, typename Dynamics, typename Measurement>
class nmhe_solve : protected nmhe_state<Scalar, NX, NU, NY, N, NC, Solver, Dynamics, Measurement>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int Ni = static_cast<int>(N);

    using base = nmhe_state<Scalar, NX, NU, NY, N, NC, Solver, Dynamics, Measurement>;

protected:
    using base::base;

    void solve_nmhe()
    {
        if(this->m_setup_failed)
        {
            fallback_to_ekf();
            return;
        }

        if(!update_arrival_cost())
        {
            fallback_to_ekf();
            return;
        }

        attempt_nmhe_solve();
    }

    void attempt_nmhe_solve()
    {
        nlp_update<Scalar> update{.x0 = this->m_warm_z};
        auto result = this->m_solver.solve(update);

        if(result.status == solve_status::optimal || result.status == solve_status::solved_inaccurate)
        {
            // A status is not a numerical postcondition; the check precedes every state and warm-start write.
            auto const primal_size = static_cast<Eigen::Index>(this->m_problem->n_vars);
            if(!detail::nlp_result_is_usable(result, result.x, primal_size))
            {
                fallback_to_ekf(solve_status::invalid_backend_result);
                return;
            }

            extract_nmhe_solution(result);
            return;
        }

        fallback_to_ekf();
    }

    bool update_arrival_cost()
    {
        auto inverse = detail::finite_full_piv_inverse(this->m_prior_cov_window[0], moving_horizon_configuration_error::non_invertible_initial_covariance);
        if(!inverse)
            return false;
        this->m_state->arrival_state = this->m_prior_state_window[0];
        this->m_state->arrival_P_inv = std::move(*inverse);
        return true;
    }

    void extract_nmhe_solution(const nlp_result<Scalar>& result)
    {
        for(std::int32_t k = 0; k <= Ni; ++k)
            this->m_x_window[static_cast<std::size_t>(k)] = result.x.segment(k * nx, nx);

        this->m_innovation = (this->m_z_window[N] - this->m_measurement(this->m_x_window[N])).eval();
        shift_warm_start(result.x);
        populate_solve_diagnostics(result);
    }

    void populate_solve_diagnostics(const nlp_result<Scalar>& result)
    {
        mhe_diagnostics<Scalar> diagnostics;
        diagnostics.status = result.status;
        diagnostics.iterations = result.iterations;
        diagnostics.solve_time = result.solve_time;
        diagnostics.cost = result.objective;
        diagnostics.primal_residual = result.primal_residual;
        diagnostics.used_ekf_fallback = false;
        this->m_diagnostics = diagnostics;

        compute_total_slack(result);
    }

    void compute_total_slack(const nlp_result<Scalar>& result)
    {
        if constexpr(NC > 0)
        {
            bool has_path_slack = this->m_config.soft_constraints && this->m_config.path_constraint_cb.has_value();
            if(has_path_slack)
            {
                Eigen::Index slack_off = (Ni + 1) * nx;
                Scalar total{0};
                for(std::int32_t k = 0; k <= Ni; ++k)
                {
                    Eigen::Map<const Vector<Scalar, NC>> sk(result.x.data() + slack_off + k * nc);
                    total += sk.sum();
                }
                this->m_diagnostics.total_slack = total;
            }
        }
    }

    void fallback_to_ekf(solve_status reason = solve_status::error)
    {
        this->m_x_window[N] = this->m_ekf.state();
        this->m_innovation = this->m_ekf.innovation();
        this->m_diagnostics = detail::ekf_fallback_diagnostics<Scalar>(reason);
    }

    // No length check: this runs only from the extraction, which the solve attempt reaches after confirming the result
    // covers the posed problem.
    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        for(std::int32_t k = 0; k < Ni; ++k)
            this->m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
        this->m_warm_z.segment(Ni * nx, nx) = this->m_ekf.state();

        shift_path_slack(sol);
    }

    void shift_path_slack(const Eigen::VectorX<Scalar>& sol)
    {
        if constexpr(NC > 0)
        {
            bool has_path_slack = this->m_config.soft_constraints && this->m_config.path_constraint_cb.has_value();
            if(has_path_slack)
            {
                Eigen::Index slack_off = (Ni + 1) * nx;
                for(std::int32_t k = 0; k < Ni; ++k)
                    this->m_warm_z.segment(slack_off + k * nc, nc) = sol.segment(slack_off + (k + 1) * nc, nc);
                this->m_warm_z.segment(slack_off + Ni * nc, nc).setZero();
            }
        }
    }
};

}

#endif
