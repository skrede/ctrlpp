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

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_qp_formulation.h"

#include "ctrlpp/ekf.h"
#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"

#include "ctrlpp/detail/numerical_diff.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"
#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/measurement_model.h"
#include "ctrlpp/mpc/differentiable_dynamics.h"
#include "ctrlpp/mpc/differentiable_measurement.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         std::size_t N,
         typename Solver,
         typename Dynamics, typename Measurement>
requires qp_solver<Solver> &&
         dynamics_model<Dynamics, Scalar, NX, NU> &&
         measurement_model<Measurement, Scalar, NX, NY>
class mhe {
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
        const mhe_config<Scalar, NX, NU, NY, N>& config)
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , ekf_{dynamics_, measurement_, ekf_config<Scalar, NX, NU, NY>{
              .Q = config.Q, .R = config.R,
              .x0 = config.x0, .P0 = config.P0,
              .numerical_eps = config.numerical_eps}}
        , arrival_cost_weight_{config.arrival_cost_weight}
        , Q_inv_{config.Q.inverse()}
        , R_inv_{config.R.inverse()}
        , eps_{config.numerical_eps}
        , x_min_{config.x_min}
        , x_max_{config.x_max}
        , residual_bound_{config.residual_bound}
        , soft_constraints_{config.soft_constraints}
        , soft_penalty_{config.soft_penalty}
        , innovation_{output_vector_t::Zero()}
    {
        // Initialize window buffers
        x_window_.fill(config.x0);
        u_window_.fill(input_vector_t::Zero());
        z_window_.fill(output_vector_t::Zero());

        // Initialize warm-start vectors
        bool has_box = x_min_.has_value() || x_max_.has_value();
        bool has_residual = residual_bound_.has_value();
        auto dims = detail::compute_mhe_dims<NX, NY>(
            N, has_box, soft_constraints_ && has_box, has_residual);

        warm_z_ = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
        warm_y_ = Eigen::VectorX<Scalar>::Zero(dims.n_con);

        // Fill initial warm-start with x0 repeated
        for (int k = 0; k <= Ni; ++k)
            warm_z_.segment(k * nx, nx) = config.x0;
    }

    void predict(const input_vector_t& u)
    {
        ekf_.predict(u);

        // Shift u_window left, insert at end
        std::rotate(u_window_.begin(), u_window_.begin() + 1, u_window_.end());
        u_window_.back() = u;

        ++step_count_;
    }

    void update(const output_vector_t& z)
    {
        ekf_.update(z);

        // Shift z_window left, insert at end
        std::rotate(z_window_.begin(), z_window_.begin() + 1, z_window_.end());
        z_window_.back() = z;

        // Warmup: state and covariance from companion EKF
        if (step_count_ < N) {
            // Update x_window from EKF during warmup
            x_window_[step_count_] = ekf_.state();
            innovation_ = ekf_.innovation();
            diagnostics_ = mhe_diagnostics<Scalar>{
                .status = solve_status::optimal,
                .used_ekf_fallback = true};
            return;
        }

        // Full MHE solve
        solve_mhe(z);
    }

    [[nodiscard]] auto state() const -> const state_vector_t& { return x_window_[N]; }
    [[nodiscard]] auto covariance() const -> const cov_matrix_t& { return ekf_.covariance(); }
    [[nodiscard]] auto innovation() const -> const output_vector_t& { return innovation_; }
    [[nodiscard]] auto trajectory() const -> std::span<const state_vector_t> { return {x_window_.data(), N + 1}; }
    [[nodiscard]] auto arrival_state() const -> const state_vector_t& { return x_window_[0]; }
    [[nodiscard]] auto arrival_covariance() const -> const cov_matrix_t& { return ekf_.covariance(); }
    [[nodiscard]] auto is_initialized() const -> bool { return step_count_ >= N; }
    [[nodiscard]] auto diagnostics() const -> const mhe_diagnostics<Scalar>& { return diagnostics_; }

private:
    void solve_mhe(const output_vector_t& z)
    {
        // Linearize dynamics around current trajectory
        auto [A_lin, B_lin] = linearize_dynamics();
        auto H_lin = linearize_measurement();

        // Get arrival cost from companion EKF
        cov_matrix_t P_arr_inv = ekf_.covariance().inverse();

        bool has_box = x_min_.has_value() || x_max_.has_value();
        bool has_residual = residual_bound_.has_value();

        // Build QP structure (P and A matrices depend on linearization)
        std::array<Matrix<Scalar, NX, NX>, 1> A_arr{A_lin};
        std::array<Matrix<Scalar, NY, NX>, 1> H_arr{H_lin};

        auto problem = detail::build_mhe_qp_structure<Scalar, NX, NU, NY>(
            N, arrival_cost_weight_, P_arr_inv, Q_inv_, R_inv_,
            A_arr, H_arr,
            has_box, soft_constraints_ && has_box, soft_penalty_,
            has_residual);

        // Build per-solve update
        std::span<const input_vector_t> u_span{u_window_.data(), N};
        std::span<const output_vector_t> z_span{z_window_.data(), N + 1};

        auto upd = detail::build_mhe_qp_update<Scalar, NX, NU, NY>(
            N, arrival_cost_weight_, P_arr_inv, Q_inv_, R_inv_,
            A_lin, B_lin, H_lin,
            ekf_.state(), u_span, z_span,
            has_box, soft_constraints_ && has_box, soft_penalty_,
            x_min_, x_max_,
            has_residual, residual_bound_,
            warm_z_, warm_y_);

        // Merge initial structure q with update q (structure has slack penalty, update has data terms)
        for (int i = 0; i < static_cast<int>(upd.q.size()); ++i)
            upd.q(i) += problem.q(i);
        problem.q = upd.q;
        problem.l = upd.l;
        problem.u = upd.u;

        // Setup and solve
        try {
            solver_.setup(problem);
            auto result = solver_.solve(
                qp_update<Scalar>{upd.q, upd.l, upd.u, warm_z_, warm_y_});

            if (result.status == solve_status::optimal ||
                result.status == solve_status::solved_inaccurate) {
                // Extract state trajectory from solution
                for (int k = 0; k <= Ni; ++k)
                    x_window_[static_cast<std::size_t>(k)] = result.x.segment(k * nx, nx);

                // Warm-start shift-and-fill for next solve
                shift_warm_start(result.x, result.y);

                // Compute innovation
                innovation_ = (z - measurement_(state())).eval();

                // Populate diagnostics
                diagnostics_ = mhe_diagnostics<Scalar>{
                    .status = result.status,
                    .iterations = result.iterations,
                    .solve_time = result.solve_time,
                    .cost = result.objective,
                    .primal_residual = result.primal_residual,
                    .dual_residual = result.dual_residual,
                    .max_constraint_violation = std::max(
                        result.primal_residual, result.dual_residual),
                    .used_ekf_fallback = false};

                // Compute slack totals if present
                auto dims = detail::compute_mhe_dims<NX, NY>(
                    N, has_box, soft_constraints_ && has_box, has_residual);
                if (dims.n_slack > 0) {
                    Scalar slack_sum{0};
                    for (int i = dims.n_states; i < dims.n_dec; ++i)
                        slack_sum += std::abs(result.x(i));
                    diagnostics_.total_slack = slack_sum;
                }

                return;
            }
        } catch (...) {
            // Solver setup or solve threw -- fall through to EKF fallback
        }

        // Fallback to companion EKF
        fallback_to_ekf();
    }

    void fallback_to_ekf()
    {
        // Fill latest window entry from EKF
        x_window_[N] = ekf_.state();
        innovation_ = ekf_.innovation();
        diagnostics_ = mhe_diagnostics<Scalar>{
            .status = solve_status::error,
            .used_ekf_fallback = true};
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol_x,
                          const Eigen::VectorX<Scalar>& sol_y)
    {
        // Shift state portion: warm[0..N-1] = sol[1..N], warm[N] = EKF prediction
        for (int k = 0; k < Ni; ++k)
            warm_z_.segment(k * nx, nx) = sol_x.segment((k + 1) * nx, nx);
        warm_z_.segment(Ni * nx, nx) = ekf_.state();

        // Shift slack portion (if present)
        auto dims = detail::compute_mhe_dims<NX, NY>(
            N, x_min_.has_value() || x_max_.has_value(),
            soft_constraints_ && (x_min_.has_value() || x_max_.has_value()),
            residual_bound_.has_value());
        if (dims.n_slack > 0) {
            int slack_off = dims.n_states;
            int per_step = nx;
            int n_steps = Ni + 1;
            for (int k = 0; k < n_steps - 1; ++k) {
                warm_z_.segment(slack_off + k * per_step, per_step) =
                    sol_x.segment(slack_off + (k + 1) * per_step, per_step);
            }
            warm_z_.segment(slack_off + (n_steps - 1) * per_step, per_step).setZero();
        }

        warm_y_ = sol_y;
    }

    [[nodiscard]] auto linearize_dynamics() const
        -> std::pair<Matrix<Scalar, NX, NX>, Matrix<Scalar, NX, NU>>
    {
        // Linearize about the midpoint of the current window trajectory
        const auto& x_ref = x_window_[N / 2];
        const auto& u_ref = u_window_[N / 2];

        Matrix<Scalar, NX, NX> A;
        Matrix<Scalar, NX, NU> B;

        if constexpr (differentiable_dynamics<Dynamics, Scalar, NX, NU>) {
            A = dynamics_.jacobian_x(x_ref, u_ref);
            B = dynamics_.jacobian_u(x_ref, u_ref);
        } else {
            A = detail::numerical_jacobian_x<Scalar, NX, NU>(dynamics_, x_ref, u_ref, eps_);
            B = detail::numerical_jacobian_u<Scalar, NX, NU>(dynamics_, x_ref, u_ref, eps_);
        }

        return {A, B};
    }

    [[nodiscard]] auto linearize_measurement() const -> Matrix<Scalar, NY, NX>
    {
        const auto& x_ref = x_window_[N / 2];

        if constexpr (differentiable_measurement<Measurement, Scalar, NX, NY>) {
            return measurement_.jacobian(x_ref);
        } else {
            return detail::numerical_jacobian_h<Scalar, NX, NY>(measurement_, x_ref, eps_);
        }
    }

    Dynamics dynamics_;
    Measurement measurement_;
    ekf<Scalar, NX, NU, NY, Dynamics, Measurement> ekf_;

    Scalar arrival_cost_weight_;
    cov_matrix_t Q_inv_;
    Matrix<Scalar, NY, NY> R_inv_;
    Scalar eps_;

    std::optional<state_vector_t> x_min_;
    std::optional<state_vector_t> x_max_;
    std::optional<output_vector_t> residual_bound_;
    bool soft_constraints_;
    Scalar soft_penalty_;

    std::array<state_vector_t, N + 1> x_window_;
    std::array<input_vector_t, N> u_window_;
    std::array<output_vector_t, N + 1> z_window_;

    std::size_t step_count_{0};
    Solver solver_;
    Eigen::VectorX<Scalar> warm_z_;
    Eigen::VectorX<Scalar> warm_y_;
    mhe_diagnostics<Scalar> diagnostics_{};
    output_vector_t innovation_;
};

// Note: No CTAD deduction guide -- Solver cannot be deduced from constructor
// arguments. Users must specify all template parameters explicitly.

namespace detail {

struct mhe_sa_dynamics {
    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>&) const -> Vector<double, 2>
    {
        return x;
    }
};

struct mhe_sa_measurement {
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        return x.template head<1>();
    }
};

struct mhe_sa_solver {
    using scalar_type = double;
    void setup(const qp_problem<double>&) {}
    auto solve(const qp_update<double>&) -> qp_result<double> { return {}; }
};

}

static_assert(ObserverPolicy<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);
static_assert(CovarianceObserver<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);

}

#endif
