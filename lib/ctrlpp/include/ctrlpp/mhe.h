#ifndef HPP_GUARD_CTRLPP_MHE_H
#define HPP_GUARD_CTRLPP_MHE_H

/// @file mhe.h
/// @brief Linear Moving Horizon Estimator with QP solver injection.
///
/// Implements constrained state estimation over a sliding window using
/// sparse QP optimization. Reuses the qp_solver concept for solver injection
/// and satisfies CovarianceObserver for drop-in replacement of EKF.
///
/// References:
///   - C.V. Rao, J.B. Rawlings, D.Q. Lee, "Constrained linear state estimation
///     -- a moving horizon approach," Automatica, 37(10), 2001.
///   - P. Kuhl, M. Diehl, T. Johansen, "Real-time optimization for large scale
///     nonlinear processes," Springer, 2011.
///   - C.V. Rao, J.B. Rawlings, J.H. Lee, "Constrained linear state estimation,"
///     Int. J. Robust Nonlinear Control, 13(10), 2003.

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/types.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_qp_formulation.h"

#include "ctrlpp/detail/numerical_diff.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"
#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"
#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"

#include <Eigen/Dense>

#include <span>
#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp {

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, typename Solver, typename Dynamics, typename Measurement>
    requires qp_solver<Solver> && dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class mhe
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr int Ni = static_cast<int>(N);

public:
    using observer_tag = struct mhe_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;

    mhe(Dynamics dynamics, Measurement measurement,
        const mhe_config<Scalar, NX, NU, NY, N> &config)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_ekf{
              m_dynamics,
              m_measurement,
              ekf_config<Scalar, NX, NU, NY>{
                  .Q             = config.Q,
                  .R             = config.R,
                  .x0            = config.x0,
                  .P0            = config.P0,
                  .numerical_eps = config.numerical_eps
              }
          }
        , m_arrival_cost_weight{config.arrival_cost_weight}
        , m_Q_inv{config.Q.inverse()}
        , m_R_inv{config.R.inverse()}
        , m_eps{config.numerical_eps}
        , m_x_min{config.x_min}
        , m_x_max{config.x_max}
        , m_residual_bound{config.residual_bound}
        , m_soft_constraints{config.soft_constraints}
        , m_soft_penalty{config.soft_penalty}
        , m_innovation{output_vector_t::Zero()}
    {
        // Initialize window buffers
        m_x_window.fill(config.x0);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());

        // Initialize warm-start vectors
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();
        auto dims = detail::compute_mhe_dims<NX, NY>(
            N, has_box, m_soft_constraints && has_box, has_residual);

        m_warm_z = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
        m_warm_y = Eigen::VectorX<Scalar>::Zero(dims.n_con);

        // Fill initial warm-start with x0 repeated
        for(int k = 0; k <= Ni; ++k)
            m_warm_z.segment(k * nx, nx) = config.x0;
    }

    void predict(const input_vector_t &u)
    {
        m_ekf.predict(u);

        // Shift u_window left, insert at end
        std::rotate(m_u_window.begin(), m_u_window.begin() + 1, m_u_window.end());
        m_u_window.back() = u;

        ++m_step_count;
    }

    void update(const output_vector_t &z)
    {
        m_ekf.update(z);

        // Shift z_window left, insert at end
        std::rotate(m_z_window.begin(), m_z_window.begin() + 1, m_z_window.end());
        m_z_window.back() = z;

        // Warmup: state and covariance from companion EKF
        if(m_step_count < N)
        {
            // Update x_window from EKF during warmup
            m_x_window[m_step_count] = m_ekf.state();
            m_innovation = m_ekf.innovation();
            m_diagnostics = mhe_diagnostics<Scalar>{
                .status            = solve_status::optimal,
                .used_ekf_fallback = true
            };
            return;
        }

        // Full MHE solve
        solve_mhe(z);
    }

    const state_vector_t &state() const
    {
        return m_x_window[N];
    }

    const cov_matrix_t &covariance() const
    {
        return m_ekf.covariance();
    }

    const output_vector_t &innovation() const
    {
        return m_innovation;
    }

    std::span<const state_vector_t> trajectory() const
    {
        return {m_x_window.data(), N + 1};
    }

    const state_vector_t &arrival_state() const
    {
        return m_x_window[0];
    }

    const cov_matrix_t &arrival_covariance() const
    {
        return m_ekf.covariance();
    }

    bool is_initialized() const
    {
        return m_step_count >= N;
    }

    const mhe_diagnostics<Scalar> &diagnostics() const
    {
        return m_diagnostics;
    }

private:
    void solve_mhe(const output_vector_t &z)
    {
        // Linearize dynamics around current trajectory
        auto [A_lin, B_lin] = linearize_dynamics();
        auto H_lin = linearize_measurement();

        // Get arrival cost from companion EKF
        cov_matrix_t P_arr_inv = m_ekf.covariance().inverse();

        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();

        // Build QP structure (P and A matrices depend on linearization)
        std::array<Matrix<Scalar, NX, NX>, 1> A_arr{A_lin};
        std::array<Matrix<Scalar, NY, NX>, 1> H_arr{H_lin};

        auto problem = detail::build_mhe_qp_structure<Scalar, NX, NU, NY>(
            N, m_arrival_cost_weight, P_arr_inv, m_Q_inv, m_R_inv,
            A_arr, H_arr,
            has_box, m_soft_constraints && has_box, m_soft_penalty,
            has_residual);

        // Build per-solve update
        std::span<const input_vector_t> u_span{m_u_window.data(), N};
        std::span<const output_vector_t> z_span{m_z_window.data(), N + 1};

        auto upd = detail::build_mhe_qp_update<Scalar, NX, NU, NY>(
            N, m_arrival_cost_weight, P_arr_inv, m_Q_inv, m_R_inv,
            A_lin, B_lin, H_lin,
            m_ekf.state(), u_span, z_span,
            has_box, m_soft_constraints && has_box, m_soft_penalty,
            m_x_min, m_x_max,
            has_residual, m_residual_bound,
            m_warm_z, m_warm_y);

        // Merge initial structure q with update q (structure has slack penalty, update has data terms)
        for(int i = 0; i < static_cast<int>(upd.q.size()); ++i)
            upd.q(i) += problem.q(i);
        problem.q = upd.q;
        problem.l = upd.l;
        problem.u = upd.u;

        // Setup and solve
        try
        {
            m_solver.setup(problem);
            auto result = m_solver.solve(
                qp_update<Scalar>{upd.q, upd.l, upd.u, m_warm_z, m_warm_y});

            if(result.status == solve_status::optimal ||
                result.status == solve_status::solved_inaccurate)
            {
                // Extract state trajectory from solution
                for(int k = 0; k <= Ni; ++k)
                    m_x_window[static_cast<std::size_t>(k)] = result.x.segment(k * nx, nx);

                // Warm-start shift-and-fill for next solve
                shift_warm_start(result.x, result.y);

                // Compute innovation
                m_innovation = (z - measurement_(state())).eval();

                // Populate diagnostics
                m_diagnostics = mhe_diagnostics<Scalar>{
                    .status                   = result.status,
                    .iterations               = result.iterations,
                    .solve_time               = result.solve_time,
                    .cost                     = result.objective,
                    .primal_residual          = result.primal_residual,
                    .dual_residual            = result.dual_residual,
                    .max_constraint_violation = std::max(
                        result.primal_residual, result.dual_residual),
                    .used_ekf_fallback = false
                };

                // Compute slack totals if present
                auto dims = detail::compute_mhe_dims<NX, NY>(
                    N, has_box, m_soft_constraints && has_box, has_residual);
                if(dims.n_slack > 0)
                {
                    Scalar slack_sum{0};
                    for(int i = dims.n_states; i < dims.n_dec; ++i)
                        slack_sum += std::abs(result.x(i));
                    m_diagnostics.total_slack = slack_sum;
                }

                return;
            }
        }
        catch(...)
        {
            // Solver setup or solve threw -- fall through to EKF fallback
        }

        // Fallback to companion EKF
        fallback_to_ekf();
    }

    void fallback_to_ekf()
    {
        // Fill latest window entry from EKF
        m_x_window[N] = m_ekf.state();
        m_innovation = m_ekf.innovation();
        m_diagnostics = mhe_diagnostics<Scalar>{
            .status            = solve_status::error,
            .used_ekf_fallback = true
        };
    }

    void shift_warm_start(const Eigen::VectorX<Scalar> &sol_x, const Eigen::VectorX<Scalar> &sol_y)
    {
        // Shift state portion: warm[0..N-1] = sol[1..N], warm[N] = EKF prediction
        for(int k = 0; k < Ni; ++k)
            m_warm_z.segment(k * nx, nx) = sol_x.segment((k + 1) * nx, nx);
        m_warm_z.segment(Ni * nx, nx) = m_ekf.state();

        // Shift slack portion (if present)
        auto dims = detail::compute_mhe_dims<NX, NY>(N, m_x_min.has_value() || m_x_max.has_value(), m_soft_constraints && (m_x_min.has_value() || m_x_max.has_value()), m_residual_bound.has_value());
        if(dims.n_slack > 0)
        {
            int slack_off = dims.n_states;
            int per_step = nx;
            int n_steps = Ni + 1;
            for(int k = 0; k < n_steps - 1; ++k)
            {
                m_warm_z.segment(slack_off + k * per_step, per_step) = sol_x.segment(slack_off + (k + 1) * per_step, per_step);
            }
            m_warm_z.segment(slack_off + (n_steps - 1) * per_step, per_step).setZero();
        }

        m_warm_y = sol_y;
    }

    std::pair<Matrix<Scalar, NX, NX>, Matrix<Scalar, NX, NU>> linearize_dynamics() const
    {
        // Linearize about the midpoint of the current window trajectory
        const auto &x_ref = m_x_window[N / 2];
        const auto &u_ref = m_u_window[N / 2];

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
        const auto &x_ref = m_x_window[N / 2];

        if constexpr(differentiable_measurement<Measurement, Scalar, NX, NY>)
            return m_measurement.jacobian(x_ref);
        else
            return detail::numerical_jacobian_h<Scalar, NX, NY>(m_measurement, x_ref, m_eps);
    }

    Dynamics m_dynamics;
    Measurement m_measurement;
    ekf<Scalar, NX, NU, NY, Dynamics, Measurement> m_ekf;

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

    std::size_t m_step_count{0};
    Solver m_solver;
    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_warm_y;
    mhe_diagnostics<Scalar> m_diagnostics{};
    output_vector_t m_innovation;
};

// Note: No CTAD deduction guide -- Solver cannot be deduced from constructor
// arguments. Users must specify all template parameters explicitly.

namespace detail {

struct mhe_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2> &x, const Vector<double, 1> &) const
    {
        return x;
    }
};

struct mhe_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2> &x) const
    {
        return x.template head<1>();
    }
};

struct mhe_sa_solver
{
    using scalar_type = double;

    void setup(const qp_problem<double> &)
    {
    }

    auto solve(const qp_update<double> &) -> qp_result<double> { return {}; }
};

}

static_assert(ObserverPolicy<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);
static_assert(CovarianceObserver<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);

}

#endif
