#ifndef HPP_GUARD_CTRLPP_NMHE_H
#define HPP_GUARD_CTRLPP_NMHE_H

/// @file nmhe.h
/// @brief Nonlinear Moving Horizon Estimator with NLP solver injection.
///
/// Implements constrained nonlinear state estimation over a sliding window
/// using multiple shooting NLP optimization. Reuses the nlp_solver concept
/// for solver injection and satisfies CovarianceObserver for drop-in
/// replacement of EKF or UKF.
///
/// References:
///   - M. Diehl, H.J. Ferreau, N. Haverbeke, "Efficient Numerical Methods
///     for Nonlinear MPC and Moving Horizon Estimation," in Nonlinear Model
///     Predictive Control, Springer, 2009.
///   - V.M. Zavala, L.T. Biegler, "The advanced-step NMPC controller:
///     Optimality, stability and robustness," Automatica, 45(1), 2009.
///   - C.V. Rao, J.B. Rawlings, D.Q. Mayne, "Constrained state estimation
///     for nonlinear discrete-time systems: stability and moving horizon
///     approximations," IEEE Trans. Automat. Control, 48(2), 2003.

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_nlp_formulation.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/types.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include <Eigen/Dense>

#include <span>
#include <array>
#include <memory>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, nlp_solver Solver, typename Dynamics, typename Measurement, std::size_t NC = 0>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class nmhe
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int Ni = static_cast<int>(N);

public:
    using observer_tag = struct nmhe_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;

    nmhe(Dynamics dynamics, Measurement measurement, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_ekf{m_dynamics, m_measurement, ekf_config<Scalar, NX, NU, NY>{.Q = config.Q, .R = config.R, .x0 = config.x0, .P0 = config.P0, .numerical_eps = config.numerical_eps}}
        , m_arrival_cost_weight{config.arrival_cost_weight}
        , m_Q_inv{config.Q.inverse()}
        , m_R_inv{config.R.inverse()}
        , m_config{config}
        , m_state{std::make_shared<nmhe_formulation_state<Scalar, NX, NU, NY, N>>()}
        , m_innovation{output_vector_t::Zero()}
    {
        initialize_buffers(config);
        build_nlp_problem();
        initialize_warm_start(config.x0);
    }

    void predict(const input_vector_t& u)
    {
        m_ekf.predict(u);
        std::rotate(m_u_window.begin(), m_u_window.begin() + 1, m_u_window.end());
        m_u_window.back() = u;
        m_state->u_window = m_u_window;
        ++m_step_count;
    }

    void update(const output_vector_t& z)
    {
        m_ekf.update(z);
        std::rotate(m_z_window.begin(), m_z_window.begin() + 1, m_z_window.end());
        m_z_window.back() = z;
        m_state->z_window = m_z_window;

        if(m_step_count < N)
        {
            m_x_window[m_step_count] = m_ekf.state();
            m_innovation = m_ekf.innovation();
            m_diagnostics = mhe_diagnostics<Scalar>{.status = solve_status::optimal, .used_ekf_fallback = true};
            return;
        }

        solve_nmhe();
    }

    const state_vector_t& state() const { return m_x_window[N]; }
    const cov_matrix_t& covariance() const { return m_ekf.covariance(); }
    const output_vector_t& innovation() const { return m_innovation; }
    std::span<const state_vector_t> trajectory() const { return {m_x_window.data(), N + 1}; }
    const state_vector_t& arrival_state() const { return m_x_window[0]; }
    const cov_matrix_t& arrival_covariance() const { return m_ekf.covariance(); }
    bool is_initialized() const { return m_step_count >= N; }
    const mhe_diagnostics<Scalar>& diagnostics() const { return m_diagnostics; }

private:
    void initialize_buffers(const nmhe_config<Scalar, NX, NU, NY, N, NC>& config)
    {
        m_x_window.fill(config.x0);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());
        m_state->arrival_state = config.x0;
        m_state->arrival_P_inv = config.P0.inverse();
    }

    void build_nlp_problem()
    {
        m_problem = detail::build_nmhe_problem<Scalar, NX, NU, NY, N, NC>(m_dynamics, m_measurement, m_config, m_state, m_Q_inv, m_R_inv);
        m_solver.setup(m_problem);
    }

    void initialize_warm_start(const state_vector_t& x0)
    {
        m_warm_z = Eigen::VectorX<Scalar>::Zero(m_problem.n_vars);
        for(int k = 0; k <= Ni; ++k)
            m_warm_z.segment(k * nx, nx) = x0;
    }

    void solve_nmhe()
    {
        update_arrival_cost();

        try
        {
            nlp_update<Scalar> update{.x0 = m_warm_z};
            auto result = m_solver.solve(update);

            if(result.status == solve_status::optimal || result.status == solve_status::solved_inaccurate)
            {
                extract_nmhe_solution(result);
                return;
            }
        }
        catch(...)
        {
        }

        fallback_to_ekf();
    }

    void update_arrival_cost()
    {
        m_state->arrival_state = m_ekf.state();
        m_state->arrival_P_inv = m_ekf.covariance().ldlt().solve(cov_matrix_t::Identity());
    }

    void extract_nmhe_solution(const nlp_result<Scalar>& result)
    {
        for(int k = 0; k <= Ni; ++k)
            m_x_window[static_cast<std::size_t>(k)] = result.x.segment(k * nx, nx);

        m_innovation = (m_z_window[N] - m_measurement(state())).eval();
        shift_warm_start(result.x);
        populate_solve_diagnostics(result);
    }

    void populate_solve_diagnostics(const nlp_result<Scalar>& result)
    {
        m_diagnostics = mhe_diagnostics<Scalar>{
            .status = result.status, .iterations = result.iterations, .solve_time = result.solve_time, .cost = result.objective, .primal_residual = result.primal_residual, .used_ekf_fallback = false};

        compute_total_slack(result);
    }

    void compute_total_slack(const nlp_result<Scalar>& result)
    {
        if constexpr(NC > 0)
        {
            bool has_path_slack = m_config.soft_constraints && m_config.path_constraint.has_value();
            if(has_path_slack)
            {
                int slack_off = (Ni + 1) * nx;
                Scalar total{0};
                for(int k = 0; k <= Ni; ++k)
                {
                    Eigen::Map<const Vector<Scalar, NC>> sk(result.x.data() + slack_off + k * nc);
                    total += sk.sum();
                }
                m_diagnostics.total_slack = total;
            }
        }
    }

    void fallback_to_ekf()
    {
        m_x_window[N] = m_ekf.state();
        m_innovation = m_ekf.innovation();
        m_diagnostics = mhe_diagnostics<Scalar>{.status = solve_status::error, .used_ekf_fallback = true};
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        for(int k = 0; k < Ni; ++k)
            m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
        m_warm_z.segment(Ni * nx, nx) = m_ekf.state();

        shift_path_slack(sol);
    }

    void shift_path_slack(const Eigen::VectorX<Scalar>& sol)
    {
        if constexpr(NC > 0)
        {
            bool has_path_slack = m_config.soft_constraints && m_config.path_constraint.has_value();
            if(has_path_slack)
            {
                int slack_off = (Ni + 1) * nx;
                for(int k = 0; k < Ni; ++k)
                    m_warm_z.segment(slack_off + k * nc, nc) = sol.segment(slack_off + (k + 1) * nc, nc);
                m_warm_z.segment(slack_off + Ni * nc, nc).setZero();
            }
        }
    }

    Dynamics m_dynamics;
    Measurement m_measurement;
    ekf<Scalar, NX, NU, NY, Dynamics, Measurement> m_ekf;

    Scalar m_arrival_cost_weight;
    cov_matrix_t m_Q_inv;
    Matrix<Scalar, NY, NY> m_R_inv;
    nmhe_config<Scalar, NX, NU, NY, N, NC> m_config;

    std::shared_ptr<nmhe_formulation_state<Scalar, NX, NU, NY, N>> m_state;
    nlp_problem<Scalar> m_problem;
    Solver m_solver{};

    std::array<state_vector_t, N + 1> m_x_window;
    std::array<input_vector_t, N> m_u_window;
    std::array<output_vector_t, N + 1> m_z_window;

    std::size_t m_step_count{0};
    Eigen::VectorX<Scalar> m_warm_z;
    mhe_diagnostics<Scalar> m_diagnostics{};
    output_vector_t m_innovation;
};

namespace detail
{

struct nmhe_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>& x, const Vector<double, 1>&) const { return x; }
};

struct nmhe_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>& x) const { return x.template head<1>(); }
};

struct nmhe_sa_solver
{
    using scalar_type = double;

    void setup(const nlp_problem<double>&) {}

    nlp_result<double> solve(const nlp_update<double>&) { return {}; }
};

}

static_assert(ObserverPolicy<nmhe<double, 2, 1, 1, 5, detail::nmhe_sa_solver, detail::nmhe_sa_dynamics, detail::nmhe_sa_measurement>>);
static_assert(CovarianceObserver<nmhe<double, 2, 1, 1, 5, detail::nmhe_sa_solver, detail::nmhe_sa_dynamics, detail::nmhe_sa_measurement>>);

}

#endif
