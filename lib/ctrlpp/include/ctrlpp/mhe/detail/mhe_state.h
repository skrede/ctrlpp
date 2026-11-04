#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_STATE_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_STATE_H

#include "ctrlpp/mhe/detail/mhe_qp_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/estimation/ekf.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"

#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"

#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <optional>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, typename Solver, typename Dynamics, typename Measurement>
class mhe_state
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int Ni = static_cast<int>(N);

protected:
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using filter_type = ekf<Scalar, NX, NU, NY, Dynamics, Measurement>;

    Dynamics m_dynamics;
    Measurement m_measurement;
    filter_type m_ekf;
    Scalar m_arrival_cost_weight;
    cov_matrix_t m_Q_inv;
    Matrix<Scalar, NY, NY> m_R_inv;
    Scalar m_eps;
    std::optional<state_vector_t> m_x_min;
    std::optional<state_vector_t> m_x_max;
    std::optional<output_vector_t> m_residual_bound;
    bool m_soft_constraints;
    Scalar m_soft_penalty;
    std::array<state_vector_t, N + 1> m_x_window;
    std::array<input_vector_t, N> m_u_window;
    std::array<output_vector_t, N + 1> m_z_window;
    std::array<state_vector_t, N + 1> m_prior_state_window;
    std::array<cov_matrix_t, N + 1> m_prior_cov_window;
    std::size_t m_step_count;
    std::size_t m_update_count;
    Solver m_solver;
    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_warm_y;
    mhe_diagnostics<Scalar> m_diagnostics;
    output_vector_t m_innovation;

    mhe_state(
            Dynamics&& dynamics, Measurement&& measurement, filter_type&& filter, const mhe_config<Scalar, NX, NU, NY, N>& config, cov_matrix_t&& q_inv, Matrix<Scalar, NY, NY>&& r_inv)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_ekf{std::move(filter)}
        , m_arrival_cost_weight{config.arrival_cost_weight}
        , m_Q_inv{std::move(q_inv)}
        , m_R_inv{std::move(r_inv)}
        , m_eps{config.numerical_eps}
        , m_x_min{config.x_min}
        , m_x_max{config.x_max}
        , m_residual_bound{config.residual_bound}
        , m_soft_constraints{config.soft_constraints}
        , m_soft_penalty{config.soft_penalty}
        , m_step_count{0}
        , m_update_count{0}
        , m_diagnostics{}
        , m_innovation{output_vector_t::Zero()}
    {
        m_x_window.fill(config.x0);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());
        m_prior_state_window.fill(config.x0);
        m_prior_cov_window.fill(config.P0);
        initialize_warm_start(config.x0);
    }

    void shift_windows(const state_vector_t& prior_state, const cov_matrix_t& prior_cov, const output_vector_t& z)
    {
        std::rotate(m_prior_state_window.begin(), m_prior_state_window.begin() + 1, m_prior_state_window.end());
        std::rotate(m_prior_cov_window.begin(), m_prior_cov_window.begin() + 1, m_prior_cov_window.end());
        m_prior_state_window.back() = prior_state;
        m_prior_cov_window.back() = prior_cov;

        std::rotate(m_z_window.begin(), m_z_window.begin() + 1, m_z_window.end());
        m_z_window.back() = z;
        ++m_update_count;
    }

    void invalidate_window_after_refusal()
    {
        auto const state = m_ekf.state();
        auto const covariance = m_ekf.covariance();
        m_x_window.fill(state);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());
        m_prior_state_window.fill(state);
        m_prior_cov_window.fill(covariance);
        m_step_count = 0;
        m_update_count = 0;
        initialize_warm_start(state);
    }

    void initialize_warm_start(const state_vector_t& x0)
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();
        auto dims = detail::compute_mhe_dims<NX, NY>(N, has_box, m_soft_constraints && has_box, has_residual);

        m_warm_z = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
        m_warm_y = Eigen::VectorX<Scalar>::Zero(dims.n_con);

        for(std::int32_t k = 0; k <= Ni; ++k)
            m_warm_z.segment(k * nx, nx) = x0;
    }

    detail::mhe_qp_dims qp_dimensions() const
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        return detail::compute_mhe_dims<NX, NY>(N, has_box, m_soft_constraints && has_box, m_residual_bound.has_value());
    }

    // A longer result is accepted: only the prefix the extraction reads and the warm start hands back must exist.
    bool result_covers_problem(const qp_result<Scalar>& result) const
    {
        auto dims = qp_dimensions();
        auto const primal_size = dims.n_dec;
        auto const dual_size = dims.n_con;
        return result.x.size() >= primal_size && result.y.size() >= dual_size && result.x.head(primal_size).allFinite() && result.y.head(dual_size).allFinite() &&
                std::isfinite(result.objective) && std::isfinite(result.solve_time) && std::isfinite(result.primal_residual) && std::isfinite(result.dual_residual) &&
                result.iterations >= 0;
    }

    void merge_structure_and_update(qp_problem<Scalar>& problem, qp_update<Scalar>& upd)
    {
        for(Eigen::Index i = 0; i < upd.q.size(); ++i)
            upd.q(i) += problem.q(i);
        problem.q = upd.q;
        problem.l = upd.l;
        problem.u = upd.u;
    }

    std::pair<Matrix<Scalar, NX, NX>, Matrix<Scalar, NX, NU>> linearize_dynamics() const
    {
        const auto& x_ref = m_x_window[N / 2];
        const auto& u_ref = m_u_window[N / 2];

        Matrix<Scalar, NX, NX> A;
        Matrix<Scalar, NX, NU> B;

        if constexpr(differentiable_dynamics<Dynamics, Scalar, NX, NU>)
        {
            A = m_dynamics.jacobian_x(x_ref, u_ref);
            B = m_dynamics.jacobian_u(x_ref, u_ref);
        }
        else
        {
            A = detail::numerical_jacobian_x<Scalar, NX, NU>(m_dynamics, x_ref, u_ref, m_eps);
            B = detail::numerical_jacobian_u<Scalar, NX, NU>(m_dynamics, x_ref, u_ref, m_eps);
        }

        return {A, B};
    }

    Matrix<Scalar, NY, NX> linearize_measurement() const
    {
        const auto& x_ref = m_x_window[N / 2];

        if constexpr(differentiable_measurement<Measurement, Scalar, NX, NY>)
            return m_measurement.jacobian(x_ref);
        else
            return detail::numerical_jacobian_h<Scalar, NX, NY>(m_measurement, x_ref, m_eps);
    }
};

}

#endif
