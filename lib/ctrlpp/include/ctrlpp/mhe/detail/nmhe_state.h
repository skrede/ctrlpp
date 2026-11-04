#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_STATE_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_STATE_H

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_solver.h"

#include "ctrlpp/estimation/ekf.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_nlp_formulation.h"

#include <Eigen/Dense>

#include <array>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC, typename Solver, typename Dynamics, typename Measurement>
class nmhe_state
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int Ni = static_cast<int>(N);

protected:
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using filter_type = ekf<Scalar, NX, NU, NY, Dynamics, Measurement>;
    using config_type = nmhe_config<Scalar, NX, NU, NY, N, NC>;

    Dynamics m_dynamics;
    Measurement m_measurement;
    filter_type m_ekf;
    Scalar m_arrival_cost_weight;
    cov_matrix_t m_Q_inv;
    Matrix<Scalar, NY, NY> m_R_inv;
    config_type m_config;
    std::shared_ptr<nmhe_formulation_state<Scalar, NX, NU, NY, N>> m_state;
    // Held at a stable heap address: a solver bridge caches the problem by pointer (nlopt_solver::setup stores &problem),
    // and the factory moves the estimator on the way out.
    std::unique_ptr<nlp_problem<Scalar>> m_problem;
    Solver m_solver;
    std::array<state_vector_t, N + 1> m_x_window;
    std::array<input_vector_t, N> m_u_window;
    std::array<output_vector_t, N + 1> m_z_window;
    std::array<state_vector_t, N + 1> m_prior_state_window;
    std::array<cov_matrix_t, N + 1> m_prior_cov_window;
    std::size_t m_step_count;
    std::size_t m_update_count;
    Eigen::VectorX<Scalar> m_warm_z;
    mhe_diagnostics<Scalar> m_diagnostics;
    output_vector_t m_innovation;
    bool m_setup_failed;

    nmhe_state(
            Dynamics&& dynamics, Measurement&& measurement, filter_type&& filter, const config_type& config, cov_matrix_t&& q_inv, Matrix<Scalar, NY, NY>&& r_inv, cov_matrix_t&& p0_inv)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_ekf{std::move(filter)}
        , m_arrival_cost_weight{config.arrival_cost_weight}
        , m_Q_inv{std::move(q_inv)}
        , m_R_inv{std::move(r_inv)}
        , m_config{config}
        , m_state{std::make_shared<nmhe_formulation_state<Scalar, NX, NU, NY, N>>()}
        , m_solver{}
        , m_step_count{0}
        , m_update_count{0}
        , m_diagnostics{}
        , m_innovation{output_vector_t::Zero()}
        , m_setup_failed{false}
    {
        initialize_buffers(config, std::move(p0_inv));
        build_nlp_problem();
        initialize_warm_start(config.x0);
    }

    void shift_priors(const state_vector_t& prior_state, const cov_matrix_t& prior_cov)
    {
        std::rotate(m_prior_state_window.begin(), m_prior_state_window.begin() + 1, m_prior_state_window.end());
        std::rotate(m_prior_cov_window.begin(), m_prior_cov_window.begin() + 1, m_prior_cov_window.end());
        m_prior_state_window.back() = prior_state;
        m_prior_cov_window.back() = prior_cov;
    }

    void shift_measurements(const output_vector_t& z)
    {
        std::rotate(m_z_window.begin(), m_z_window.begin() + 1, m_z_window.end());
        m_z_window.back() = z;
        m_state->z_window = m_z_window;
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
        m_state->arrival_state = state;
        m_state->u_window = m_u_window;
        m_state->z_window = m_z_window;
        m_step_count = 0;
        m_update_count = 0;
        initialize_warm_start(state);
    }

    void initialize_buffers(const config_type& config, cov_matrix_t p0_inv)
    {
        m_x_window.fill(config.x0);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());
        m_prior_state_window.fill(config.x0);
        m_prior_cov_window.fill(config.P0);
        m_state->arrival_state = config.x0;
        m_state->arrival_P_inv = std::move(p0_inv);
    }

    void build_nlp_problem()
    {
        m_problem = std::make_unique<nlp_problem<Scalar>>(detail::build_nmhe_problem<Scalar, NX, NU, NY, N, NC>(m_dynamics, m_measurement, m_config, m_state, m_Q_inv, m_R_inv));
        m_setup_failed = !detail::setup_nlp_solver(m_solver, *m_problem).has_value();
    }

    void initialize_warm_start(const state_vector_t& x0)
    {
        m_warm_z = Eigen::VectorX<Scalar>::Zero(m_problem->n_vars);
        for(std::int32_t k = 0; k <= Ni; ++k)
            m_warm_z.segment(k * nx, nx) = x0;
    }
};

}

#endif
