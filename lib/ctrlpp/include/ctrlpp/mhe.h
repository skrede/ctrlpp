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

namespace ctrlpp
{

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

    mhe(Dynamics dynamics, Measurement measurement, const mhe_config<Scalar, NX, NU, NY, N>& config)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_ekf{m_dynamics, m_measurement, ekf_config<Scalar, NX, NU, NY>{.Q = config.Q, .R = config.R, .x0 = config.x0, .P0 = config.P0, .numerical_eps = config.numerical_eps}}
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
        m_x_window.fill(config.x0);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());
        initialize_warm_start(config.x0);
    }

    void predict(const input_vector_t& u)
    {
        m_ekf.predict(u);
        std::rotate(m_u_window.begin(), m_u_window.begin() + 1, m_u_window.end());
        m_u_window.back() = u;
        ++m_step_count;
    }

    void update(const output_vector_t& z)
    {
        m_ekf.update(z);
        std::rotate(m_z_window.begin(), m_z_window.begin() + 1, m_z_window.end());
        m_z_window.back() = z;

        if(m_step_count < N)
        {
            m_x_window[m_step_count] = m_ekf.state();
            m_innovation = m_ekf.innovation();
            m_diagnostics = mhe_diagnostics<Scalar>{.status = solve_status::optimal, .used_ekf_fallback = true};
            return;
        }

        solve_mhe(z);
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
    void initialize_warm_start(const state_vector_t& x0)
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();
        auto dims = detail::compute_mhe_dims<NX, NY>(N, has_box, m_soft_constraints && has_box, has_residual);

        m_warm_z = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
        m_warm_y = Eigen::VectorX<Scalar>::Zero(dims.n_con);

        for(int k = 0; k <= Ni; ++k)
            m_warm_z.segment(k * nx, nx) = x0;
    }

    void solve_mhe(const output_vector_t& z)
    {
        auto [A_lin, B_lin] = linearize_dynamics();
        auto H_lin = linearize_measurement();
        auto problem = build_qp_structure(A_lin, H_lin);
        auto upd = build_qp_update(A_lin, B_lin, H_lin);

        merge_structure_and_update(problem, upd);

        try
        {
            m_solver.setup(problem);
            auto result = m_solver.solve(qp_update<Scalar>{upd.q, upd.l, upd.u, m_warm_z, m_warm_y});

            if(result.status == solve_status::optimal || result.status == solve_status::solved_inaccurate)
            {
                extract_mhe_solution(result, z);
                return;
            }
        }
        catch(...)
        {
        }

        fallback_to_ekf();
    }

    auto build_qp_structure(const Matrix<Scalar, NX, NX>& A_lin, const Matrix<Scalar, NY, NX>& H_lin) -> qp_problem<Scalar>
    {
        cov_matrix_t P_arr_inv = m_ekf.covariance().inverse();
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();
        std::array<Matrix<Scalar, NX, NX>, 1> A_arr{A_lin};
        std::array<Matrix<Scalar, NY, NX>, 1> H_arr{H_lin};

        return detail::build_mhe_qp_structure<Scalar, NX, NU, NY>(N, m_arrival_cost_weight, P_arr_inv, m_Q_inv, m_R_inv, A_arr, H_arr, has_box, m_soft_constraints && has_box, m_soft_penalty, has_residual);
    }

    auto build_qp_update(const Matrix<Scalar, NX, NX>& A_lin, const Matrix<Scalar, NX, NU>& B_lin, const Matrix<Scalar, NY, NX>& H_lin) -> qp_update<Scalar>
    {
        cov_matrix_t P_arr_inv = m_ekf.covariance().inverse();
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();
        std::span<const input_vector_t> u_span{m_u_window.data(), N};
        std::span<const output_vector_t> z_span{m_z_window.data(), N + 1};

        return detail::build_mhe_qp_update<Scalar, NX, NU, NY>(N, m_arrival_cost_weight, P_arr_inv, m_Q_inv, m_R_inv, A_lin, B_lin, H_lin, m_ekf.state(), u_span, z_span, has_box, m_soft_constraints && has_box, m_soft_penalty, m_x_min, m_x_max, has_residual, m_residual_bound, m_warm_z, m_warm_y);
    }

    void merge_structure_and_update(qp_problem<Scalar>& problem, qp_update<Scalar>& upd)
    {
        for(int i = 0; i < static_cast<int>(upd.q.size()); ++i)
            upd.q(i) += problem.q(i);
        problem.q = upd.q;
        problem.l = upd.l;
        problem.u = upd.u;
    }

    void extract_mhe_solution(const qp_result<Scalar>& result, const output_vector_t& z)
    {
        for(int k = 0; k <= Ni; ++k)
            m_x_window[static_cast<std::size_t>(k)] = result.x.segment(k * nx, nx);

        shift_warm_start(result.x, result.y);
        m_innovation = (z - m_measurement(state())).eval();
        populate_solve_diagnostics(result);
    }

    void populate_solve_diagnostics(const qp_result<Scalar>& result)
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();

        m_diagnostics = mhe_diagnostics<Scalar>{.status = result.status,
                                                .iterations = result.iterations,
                                                .solve_time = result.solve_time,
                                                .cost = result.objective,
                                                .primal_residual = result.primal_residual,
                                                .dual_residual = result.dual_residual,
                                                .max_constraint_violation = std::max(result.primal_residual, result.dual_residual),
                                                .used_ekf_fallback = false};

        auto dims = detail::compute_mhe_dims<NX, NY>(N, has_box, m_soft_constraints && has_box, has_residual);
        if(dims.n_slack > 0)
        {
            Scalar slack_sum{0};
            for(int i = dims.n_states; i < dims.n_dec; ++i)
                slack_sum += std::abs(result.x(i));
            m_diagnostics.total_slack = slack_sum;
        }
    }

    void fallback_to_ekf()
    {
        m_x_window[N] = m_ekf.state();
        m_innovation = m_ekf.innovation();
        m_diagnostics = mhe_diagnostics<Scalar>{.status = solve_status::error, .used_ekf_fallback = true};
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol_x, const Eigen::VectorX<Scalar>& sol_y)
    {
        for(int k = 0; k < Ni; ++k)
            m_warm_z.segment(k * nx, nx) = sol_x.segment((k + 1) * nx, nx);
        m_warm_z.segment(Ni * nx, nx) = m_ekf.state();

        shift_slack_warm_start(sol_x);
        m_warm_y = sol_y;
    }

    void shift_slack_warm_start(const Eigen::VectorX<Scalar>& sol_x)
    {
        auto dims = detail::compute_mhe_dims<NX, NY>(N, m_x_min.has_value() || m_x_max.has_value(), m_soft_constraints && (m_x_min.has_value() || m_x_max.has_value()), m_residual_bound.has_value());
        if(dims.n_slack <= 0)
            return;

        int slack_off = dims.n_states;
        int per_step = nx;
        int n_steps = Ni + 1;
        for(int k = 0; k < n_steps - 1; ++k)
            m_warm_z.segment(slack_off + k * per_step, per_step) = sol_x.segment(slack_off + (k + 1) * per_step, per_step);
        m_warm_z.segment(slack_off + (n_steps - 1) * per_step, per_step).setZero();
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

namespace detail
{

struct mhe_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>& x, const Vector<double, 1>&) const { return x; }
};

struct mhe_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>& x) const { return x.template head<1>(); }
};

struct mhe_sa_solver
{
    using scalar_type = double;

    void setup(const qp_problem<double>&) {}

    auto solve(const qp_update<double>&) -> qp_result<double> { return {}; }
};

}

static_assert(ObserverPolicy<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);
static_assert(CovarianceObserver<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);

}

#endif
