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

#include "ctrlpp/ekf.h"
#include "ctrlpp/observer_policy.h"
#include "ctrlpp/types.h"

#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/measurement_model.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <utility>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         std::size_t N,
         nlp_solver Solver,
         typename Dynamics, typename Measurement,
         std::size_t NC = 0>
requires dynamics_model<Dynamics, Scalar, NX, NU> &&
         measurement_model<Measurement, Scalar, NX, NY>
class nmhe {
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

    nmhe(Dynamics dynamics, Measurement measurement,
         const nmhe_config<Scalar, NX, NU, NY, N, NC>& config)
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , ekf_{dynamics_, measurement_, ekf_config<Scalar, NX, NU, NY>{
              .Q = config.Q, .R = config.R,
              .x0 = config.x0, .P0 = config.P0,
              .numerical_eps = config.numerical_eps}}
        , arrival_cost_weight_{config.arrival_cost_weight}
        , Q_inv_{config.Q.inverse()}
        , R_inv_{config.R.inverse()}
        , config_{config}
        , state_{std::make_shared<nmhe_formulation_state<Scalar, NX, NU, NY, N>>()}
        , innovation_{output_vector_t::Zero()}
    {
        // Initialize window buffers
        x_window_.fill(config.x0);
        u_window_.fill(input_vector_t::Zero());
        z_window_.fill(output_vector_t::Zero());

        // Initialize formulation state
        state_->arrival_state = config.x0;
        state_->arrival_P_inv = config.P0.inverse();

        // Build NLP problem
        problem_ = detail::build_nmhe_problem<Scalar, NX, NU, NY, N, NC>(
            dynamics_, measurement_, config_, state_, Q_inv_, R_inv_);

        solver_.setup(problem_);

        // Pre-allocate warm-start vector
        warm_z_ = Eigen::VectorX<Scalar>::Zero(problem_.n_vars);
        for (int k = 0; k <= Ni; ++k) {
            warm_z_.segment(k * nx, nx) = config.x0;
        }
    }

    void predict(const input_vector_t& u)
    {
        ekf_.predict(u);

        // Shift u_window left, insert at end
        std::rotate(u_window_.begin(), u_window_.begin() + 1, u_window_.end());
        u_window_.back() = u;

        // Mirror to formulation state
        state_->u_window = u_window_;

        ++step_count_;
    }

    void update(const output_vector_t& z)
    {
        ekf_.update(z);

        // Shift z_window left, insert at end
        std::rotate(z_window_.begin(), z_window_.begin() + 1, z_window_.end());
        z_window_.back() = z;

        // Mirror to formulation state
        state_->z_window = z_window_;

        // Warmup: use companion EKF
        if (step_count_ < N) {
            x_window_[step_count_] = ekf_.state();
            innovation_ = ekf_.innovation();
            diagnostics_ = mhe_diagnostics<Scalar>{
                .status = solve_status::optimal,
                .used_ekf_fallback = true};
            return;
        }

        // Full MHE NLP solve
        solve_nmhe();
    }

    [[nodiscard]] auto state() const -> const state_vector_t& { return x_window_[N]; }
    [[nodiscard]] auto covariance() const -> const cov_matrix_t& { return ekf_.covariance(); }
    [[nodiscard]] auto innovation() const -> const output_vector_t& { return innovation_; }

    [[nodiscard]] auto trajectory() const -> std::span<const state_vector_t>
    {
        return {x_window_.data(), N + 1};
    }

    [[nodiscard]] auto arrival_state() const -> const state_vector_t& { return x_window_[0]; }
    [[nodiscard]] auto arrival_covariance() const -> const cov_matrix_t& { return ekf_.covariance(); }
    [[nodiscard]] auto is_initialized() const -> bool { return step_count_ >= N; }
    [[nodiscard]] auto diagnostics() const -> const mhe_diagnostics<Scalar>& { return diagnostics_; }

private:
    void solve_nmhe()
    {
        // Update arrival cost from companion EKF
        state_->arrival_state = ekf_.state();
        state_->arrival_P_inv = ekf_.covariance().ldlt().solve(
            cov_matrix_t::Identity());

        // Warm-start solve
        nlp_update<Scalar> update{.x0 = warm_z_};

        try {
            auto result = solver_.solve(update);

            if (result.status == solve_status::optimal ||
                result.status == solve_status::solved_inaccurate) {
                // Extract state trajectory from solution
                for (int k = 0; k <= Ni; ++k) {
                    x_window_[static_cast<std::size_t>(k)] =
                        result.x.segment(k * nx, nx);
                }

                // Compute innovation from latest estimate
                innovation_ = (z_window_[N] - measurement_(state())).eval();

                // Warm-start shift-and-fill for next solve
                shift_warm_start(result.x);

                // Populate diagnostics
                diagnostics_ = mhe_diagnostics<Scalar>{
                    .status = result.status,
                    .iterations = result.iterations,
                    .solve_time = result.solve_time,
                    .cost = result.objective,
                    .primal_residual = result.primal_residual,
                    .used_ekf_fallback = false};

                // Compute total slack
                if constexpr (NC > 0) {
                    bool has_path_slack = config_.soft_constraints
                        && config_.path_constraint.has_value();
                    if (has_path_slack) {
                        int slack_off = (Ni + 1) * nx;
                        Scalar total{0};
                        for (int k = 0; k <= Ni; ++k) {
                            Eigen::Map<const Vector<Scalar, NC>> sk(
                                result.x.data() + slack_off + k * nc);
                            total += sk.sum();
                        }
                        diagnostics_.total_slack = total;
                    }
                }

                return;
            }
        } catch (...) {
            // Solver threw -- fall through to EKF fallback
        }

        // Fallback to companion EKF
        fallback_to_ekf();
    }

    void fallback_to_ekf()
    {
        x_window_[N] = ekf_.state();
        innovation_ = ekf_.innovation();
        diagnostics_ = mhe_diagnostics<Scalar>{
            .status = solve_status::error,
            .used_ekf_fallback = true};
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        // Shift states: warm[k] = sol[k+1] for k=0..N-1, warm[N] = EKF prediction
        for (int k = 0; k < Ni; ++k) {
            warm_z_.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
        }
        warm_z_.segment(Ni * nx, nx) = ekf_.state();

        // Shift path slack (if present)
        if constexpr (NC > 0) {
            bool has_path_slack = config_.soft_constraints
                && config_.path_constraint.has_value();
            if (has_path_slack) {
                int slack_off = (Ni + 1) * nx;
                for (int k = 0; k < Ni; ++k) {
                    warm_z_.segment(slack_off + k * nc, nc) =
                        sol.segment(slack_off + (k + 1) * nc, nc);
                }
                warm_z_.segment(slack_off + Ni * nc, nc).setZero();
            }
        }
    }

    Dynamics dynamics_;
    Measurement measurement_;
    ekf<Scalar, NX, NU, NY, Dynamics, Measurement> ekf_;

    Scalar arrival_cost_weight_;
    cov_matrix_t Q_inv_;
    Matrix<Scalar, NY, NY> R_inv_;
    nmhe_config<Scalar, NX, NU, NY, N, NC> config_;

    std::shared_ptr<nmhe_formulation_state<Scalar, NX, NU, NY, N>> state_;
    nlp_problem<Scalar> problem_;
    Solver solver_{};

    std::array<state_vector_t, N + 1> x_window_;
    std::array<input_vector_t, N> u_window_;
    std::array<output_vector_t, N + 1> z_window_;

    std::size_t step_count_{0};
    Eigen::VectorX<Scalar> warm_z_;
    mhe_diagnostics<Scalar> diagnostics_{};
    output_vector_t innovation_;
};

// Note: No CTAD deduction guide -- Solver cannot be deduced from constructor
// arguments. Users must specify all template parameters explicitly.

namespace detail {

struct nmhe_sa_dynamics {
    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>&) const -> Vector<double, 2>
    {
        return x;
    }
};

struct nmhe_sa_measurement {
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        return x.template head<1>();
    }
};

struct nmhe_sa_solver {
    using scalar_type = double;
    void setup(const nlp_problem<double>&) {}
    auto solve(const nlp_update<double>&) -> nlp_result<double> { return {}; }
};

}

static_assert(ObserverPolicy<nmhe<double, 2, 1, 1, 5, detail::nmhe_sa_solver, detail::nmhe_sa_dynamics, detail::nmhe_sa_measurement>>);
static_assert(CovarianceObserver<nmhe<double, 2, 1, 1, 5, detail::nmhe_sa_solver, detail::nmhe_sa_dynamics, detail::nmhe_sa_measurement>>);

}

#endif
