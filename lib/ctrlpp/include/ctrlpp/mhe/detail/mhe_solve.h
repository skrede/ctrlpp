#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_SOLVE_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_SOLVE_H

#include "ctrlpp/mhe/detail/mhe_state.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_qp_formulation.h"

#include <Eigen/Dense>

#include <span>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, typename Solver, typename Dynamics, typename Measurement>
class mhe_solve : protected mhe_state<Scalar, NX, NU, NY, N, Solver, Dynamics, Measurement>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int Ni = static_cast<int>(N);

    using base = mhe_state<Scalar, NX, NU, NY, N, Solver, Dynamics, Measurement>;

protected:
    using typename base::input_vector_t;
    using typename base::output_vector_t;
    using typename base::cov_matrix_t;

    using base::base;

    void use_filter_estimate()
    {
        if(this->m_step_count < N)
            this->m_x_window[this->m_step_count] = this->m_ekf.state();
        this->m_x_window[N] = this->m_ekf.state();
        this->m_innovation = this->m_ekf.innovation();
        this->m_diagnostics = detail::ekf_fallback_diagnostics<Scalar>(solve_status::optimal);
    }

    void solve_mhe(const output_vector_t& z)
    {
        auto arrival_inverse = detail::finite_full_piv_inverse(this->m_prior_cov_window[0], moving_horizon_configuration_error::non_invertible_initial_covariance);
        if(!arrival_inverse)
        {
            fallback_to_ekf();
            return;
        }

        auto [A_lin, B_lin] = this->linearize_dynamics();
        auto H_lin = this->linearize_measurement();
        auto problem = build_qp_structure(A_lin, H_lin, *arrival_inverse);
        auto upd = build_qp_update(A_lin, B_lin, H_lin, *arrival_inverse);

        this->merge_structure_and_update(problem, upd);

        attempt_mhe_solve(problem, upd, z);
    }

    void attempt_mhe_solve(const qp_problem<Scalar>& problem, const qp_update<Scalar>& upd, const output_vector_t& z)
    {
        if(!detail::setup_qp_solver(this->m_solver, problem).has_value())
        {
            fallback_to_ekf();
            return;
        }

        auto result = this->m_solver.solve(qp_update<Scalar>{upd.q, upd.l, upd.u, this->m_warm_z, this->m_warm_y});

        if(result.status == solve_status::optimal || result.status == solve_status::solved_inaccurate)
        {
            // A status is not a shape: an accepted result may still be too short for the window writes that follow.
            if(!this->result_covers_problem(result))
            {
                fallback_to_ekf(solve_status::invalid_backend_result);
                return;
            }

            extract_mhe_solution(result, z);
            return;
        }

        fallback_to_ekf();
    }

    qp_problem<Scalar> build_qp_structure(const Matrix<Scalar, NX, NX>& A_lin, const Matrix<Scalar, NY, NX>& H_lin, const cov_matrix_t& arrival_inverse)
    {
        bool has_box = this->m_x_min.has_value() || this->m_x_max.has_value();
        bool has_residual = this->m_residual_bound.has_value();
        std::array<Matrix<Scalar, NX, NX>, 1> A_arr{A_lin};
        std::array<Matrix<Scalar, NY, NX>, 1> H_arr{H_lin};

        return detail::build_mhe_qp_structure<Scalar, NX, NU, NY>(N, this->m_arrival_cost_weight, arrival_inverse, this->m_Q_inv, this->m_R_inv, A_arr, H_arr, has_box,
                this->m_soft_constraints && has_box, this->m_soft_penalty, has_residual);
    }

    qp_update<Scalar> build_qp_update(const Matrix<Scalar, NX, NX>& A_lin, const Matrix<Scalar, NX, NU>& B_lin, const Matrix<Scalar, NY, NX>& H_lin, const cov_matrix_t& arrival_inverse)
    {
        bool has_box = this->m_x_min.has_value() || this->m_x_max.has_value();
        std::span<const input_vector_t> u_span{this->m_u_window.data(), N};
        std::span<const output_vector_t> z_span{this->m_z_window.data(), N + 1};

        return detail::build_mhe_qp_update<Scalar, NX, NU, NY>(N, this->m_arrival_cost_weight, arrival_inverse, this->m_Q_inv, this->m_R_inv, A_lin, B_lin, H_lin,
                this->m_prior_state_window[0], u_span, z_span, has_box, this->m_soft_constraints && has_box, this->m_soft_penalty, this->m_x_min, this->m_x_max, this->m_residual_bound,
                this->m_warm_z, this->m_warm_y);
    }

    void extract_mhe_solution(const qp_result<Scalar>& result, const output_vector_t& z)
    {
        for(std::int32_t k = 0; k <= Ni; ++k)
            this->m_x_window[static_cast<std::size_t>(k)] = result.x.segment(k * nx, nx);

        shift_warm_start(result.x, result.y);
        this->m_innovation = (z - this->m_measurement(this->m_x_window[N])).eval();
        populate_solve_diagnostics(result);
    }

    static mhe_diagnostics<Scalar> solve_diagnostics(const qp_result<Scalar>& result)
    {
        mhe_diagnostics<Scalar> diagnostics;
        diagnostics.status = result.status;
        diagnostics.iterations = result.iterations;
        diagnostics.solve_time = result.solve_time;
        diagnostics.cost = result.objective;
        diagnostics.primal_residual = result.primal_residual;
        diagnostics.dual_residual = result.dual_residual;
        diagnostics.max_constraint_violation = std::max(result.primal_residual, result.dual_residual);
        diagnostics.used_ekf_fallback = false;
        return diagnostics;
    }

    void populate_solve_diagnostics(const qp_result<Scalar>& result)
    {
        bool has_box = this->m_x_min.has_value() || this->m_x_max.has_value();
        bool has_residual = this->m_residual_bound.has_value();

        this->m_diagnostics = solve_diagnostics(result);

        auto dims = detail::compute_mhe_dims<NX, NY>(N, has_box, this->m_soft_constraints && has_box, has_residual);
        if(dims.n_slack > 0)
        {
            Scalar slack_sum{0};
            for(Eigen::Index i = dims.n_states; i < dims.n_dec; ++i)
                slack_sum += std::abs(result.x(i));
            this->m_diagnostics.total_slack = slack_sum;
        }
    }

    void fallback_to_ekf(solve_status reason = solve_status::error)
    {
        this->m_x_window[N] = this->m_ekf.state();
        this->m_innovation = this->m_ekf.innovation();
        this->m_diagnostics = detail::ekf_fallback_diagnostics<Scalar>(reason);
    }

    // No length check of its own: only the extraction calls this, after the result was found to cover the problem.
    void shift_warm_start(const Eigen::VectorX<Scalar>& sol_x, const Eigen::VectorX<Scalar>& sol_y)
    {
        for(std::int32_t k = 0; k < Ni; ++k)
            this->m_warm_z.segment(k * nx, nx) = sol_x.segment((k + 1) * nx, nx);
        this->m_warm_z.segment(Ni * nx, nx) = this->m_ekf.state();

        shift_slack_warm_start(sol_x);
        this->m_warm_y = sol_y;
    }

    void shift_slack_warm_start(const Eigen::VectorX<Scalar>& /*sol_x*/)
    {
        auto dims = detail::compute_mhe_dims<NX, NY>(N, this->m_x_min.has_value() || this->m_x_max.has_value(),
                this->m_soft_constraints && (this->m_x_min.has_value() || this->m_x_max.has_value()), this->m_residual_bound.has_value());
        if(dims.n_slack <= 0)
            return;

        // Non-negative slacks warm-start at zero: the softened bound is inactive unless the shifted estimate reenters
        // the constraint region.
        this->m_warm_z.segment(dims.n_states, dims.n_slack).setZero();
    }
};

}

#endif
