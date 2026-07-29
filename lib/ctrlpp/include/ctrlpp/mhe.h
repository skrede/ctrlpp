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

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_qp_formulation.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"
#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"

#include "ctrlpp/detail/numerical_diff.h"

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
    static_assert(N > 0, "Window length N must be positive: it sizes the fixed estimation window arrays, which the update rotates, reads the trailing element of, and indexes at their midpoint, none of which is defined for an empty window");

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

    /// @brief Fallible factory, and the only way to originate an estimator.
    ///
    /// The error is a variant: embedded-filter configuration failures retain
    /// their `filter_error` type, while formulation-specific failures use
    /// `moving_horizon_configuration_error`. The latter rejects singular Q, R,
    /// or P0 before solving for inverse weights, and validates every constraint,
    /// penalty, and finite-difference operand consumed by the QP.
    static auto create(Dynamics dynamics, Measurement measurement,
                       const mhe_config<Scalar, NX, NU, NY, N>& config)
        -> ctrlpp::expected<mhe, moving_horizon_construction_error>
    {
        auto filter = ekf<Scalar, NX, NU, NY, Dynamics, Measurement>::create(
            dynamics, measurement, ekf_config<Scalar, NX, NU, NY>{.Q = config.Q, .R = config.R, .x0 = config.x0, .P0 = config.P0, .numerical_eps = config.numerical_eps});
        if(!filter)
            return ctrlpp::unexpected(
                moving_horizon_construction_error{filter.error()});

        auto q_inv = detail::finite_full_piv_inverse(
            config.Q,
            moving_horizon_configuration_error::
                non_invertible_process_noise);
        if(!q_inv)
            return ctrlpp::unexpected(
                moving_horizon_construction_error{q_inv.error()});
        auto r_inv = detail::finite_full_piv_inverse(
            config.R,
            moving_horizon_configuration_error::
                non_invertible_measurement_noise);
        if(!r_inv)
            return ctrlpp::unexpected(
                moving_horizon_construction_error{r_inv.error()});
        auto p0_inv = detail::finite_full_piv_inverse(
            config.P0,
            moving_horizon_configuration_error::
                non_invertible_initial_covariance);
        if(!p0_inv)
            return ctrlpp::unexpected(moving_horizon_construction_error{
                p0_inv.error()});
        auto options = detail::validate_moving_horizon_options(config);
        if(!options)
            return ctrlpp::unexpected(
                moving_horizon_construction_error{options.error()});

        return mhe{validated_tag{},
                   std::move(dynamics),
                   std::move(measurement),
                   std::move(*filter),
                   config,
                   std::move(*q_inv),
                   std::move(*r_inv)};
    }

private:
    /// @brief Tag selecting the non-validating constructor reserved for
    /// `create`, which is what makes the factory the only public path and the
    /// validation impossible to bypass.
    struct validated_tag
    {
    };

    mhe(validated_tag, Dynamics dynamics, Measurement measurement,
        ekf<Scalar, NX, NU, NY, Dynamics, Measurement> filter,
        const mhe_config<Scalar, NX, NU, NY, N>& config,
        cov_matrix_t q_inv,
        Matrix<Scalar, NY, NY> r_inv)
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
        , m_innovation{output_vector_t::Zero()}
    {
        m_x_window.fill(config.x0);
        m_u_window.fill(input_vector_t::Zero());
        m_z_window.fill(output_vector_t::Zero());
        m_prior_state_window.fill(config.x0);
        m_prior_cov_window.fill(config.P0);
        initialize_warm_start(config.x0);
    }

public:
    void predict(const input_vector_t& u)
    {
        m_ekf.predict(u);
        std::rotate(m_u_window.begin(), m_u_window.begin() + 1, m_u_window.end());
        m_u_window.back() = u;
        ++m_step_count;
    }

    /// @brief Incorporate a measurement, or report why it could not be.
    ///
    /// The embedded filter classifies the measurement against the carried
    /// estimate, and this type FORWARDS that verdict verbatim rather than naming
    /// the same three conditions a second time. The estimator has no rejection of
    /// its own -- every operand it could refuse is one the filter already
    /// refuses -- and a second enumeration would be free to drift from the one
    /// that actually decides.
    ///
    /// A refusal cannot undo the prediction: its input was already applied to
    /// the plant, so the embedded filter keeps that finite prior. The fixed-step
    /// window cannot represent the missing measurement, however, and must not
    /// retain the input beside an older measurement history. The estimator
    /// therefore invalidates the complete window and uses the embedded filter
    /// until N new accepted measurements refill a coherent horizon.
    ///
    /// `diagnostics()` stays what it is: a report on a step that SUCCEEDED --
    /// which of the two estimators produced the estimate, and how the solve
    /// went. After a refusal it still describes the last accepted measurement,
    /// while `state()` exposes the current predicted prior. The return value is
    /// what distinguishes that uncorrected prior from an accepted estimate.
    auto update(const output_vector_t& z) -> ctrlpp::expected<void, ekf_update_error>
    {
        // Capture the predicted (prior) estimate at the current time, before the
        // measurement correction. The window shift below places it so the oldest
        // entry holds the estimate at the window head formed only from data
        // strictly before the window, which is the arrival prior that avoids
        // double-counting the window measurements.
        const state_vector_t prior_state = m_ekf.state();
        const cov_matrix_t prior_cov = m_ekf.covariance();

        if(const auto stepped = m_ekf.update(z); !stepped)
        {
            invalidate_window_after_refusal();
            return ctrlpp::unexpected(stepped.error());
        }

        std::rotate(m_prior_state_window.begin(), m_prior_state_window.begin() + 1, m_prior_state_window.end());
        std::rotate(m_prior_cov_window.begin(), m_prior_cov_window.begin() + 1, m_prior_cov_window.end());
        m_prior_state_window.back() = prior_state;
        m_prior_cov_window.back() = prior_cov;

        std::rotate(m_z_window.begin(), m_z_window.begin() + 1, m_z_window.end());
        m_z_window.back() = z;
        ++m_update_count;

        if(m_update_count <= N)
        {
            if(m_step_count < N)
                m_x_window[m_step_count] = m_ekf.state();
            m_x_window[N] = m_ekf.state();
            m_innovation = m_ekf.innovation();
            m_diagnostics = mhe_diagnostics<Scalar>{.status = solve_status::optimal, .used_ekf_fallback = true};
            return {};
        }

        solve_mhe(z);
        return {};
    }

    const state_vector_t& state() const { return m_x_window[N]; }
    const cov_matrix_t& covariance() const { return m_ekf.covariance(); }
    const output_vector_t& innovation() const { return m_innovation; }
    std::span<const state_vector_t> trajectory() const { return {m_x_window.data(), N + 1}; }
    const state_vector_t& arrival_state() const { return m_x_window[0]; }
    const cov_matrix_t& arrival_covariance() const { return m_prior_cov_window[0]; }
    bool is_initialized() const { return m_update_count > N; }
    const mhe_diagnostics<Scalar>& diagnostics() const { return m_diagnostics; }

private:
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

        for(int k = 0; k <= Ni; ++k)
            m_warm_z.segment(k * nx, nx) = x0;
    }

    void solve_mhe(const output_vector_t& z)
    {
        auto arrival_inverse =
            detail::finite_full_piv_inverse(
                m_prior_cov_window[0],
                moving_horizon_configuration_error::
                    non_invertible_initial_covariance);
        if(!arrival_inverse)
        {
            fallback_to_ekf();
            return;
        }

        auto [A_lin, B_lin] = linearize_dynamics();
        auto H_lin = linearize_measurement();
        auto problem = build_qp_structure(A_lin, H_lin, *arrival_inverse);
        auto upd =
            build_qp_update(A_lin, B_lin, H_lin, *arrival_inverse);

        merge_structure_and_update(problem, upd);

        attempt_mhe_solve(problem, upd, z);
    }

    void attempt_mhe_solve(const qp_problem<Scalar>& problem, const qp_update<Scalar>& upd, const output_vector_t& z)
    {
        if(!detail::setup_qp_solver(m_solver, problem).has_value())
        {
            fallback_to_ekf();
            return;
        }

        auto result = m_solver.solve(qp_update<Scalar>{upd.q, upd.l, upd.u, m_warm_z, m_warm_y});

        if(result.status == solve_status::optimal || result.status == solve_status::solved_inaccurate)
        {
            // A status is not a shape: the accept-set above is decided purely
            // from what the backend reports, so an accepted result may still be
            // too short for the window writes that follow. Checked here, before
            // the extraction, and reported through the same channel a setup
            // failure uses.
            if(!result_covers_problem(result))
            {
                fallback_to_ekf(solve_status::invalid_backend_result);
                return;
            }

            extract_mhe_solution(result, z);
            return;
        }

        fallback_to_ekf();
    }

    /// @brief Dimensions of the QP this estimator poses, derived from the window
    /// length and whichever optional bound blocks the configuration carries.
    auto qp_dimensions() const -> detail::mhe_qp_dims
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        return detail::compute_mhe_dims<NX, NY>(N, has_box, m_soft_constraints && has_box, m_residual_bound.has_value());
    }

    /// @brief Whether a backend result is complete and finite for everything
    /// the extraction and diagnostics consume.
    ///
    /// The primal is sliced per window node at offsets derived from the window
    /// length, and the dual is stored and handed straight back to the backend as
    /// the next warm start, so both are compared against the dimensions this
    /// estimator derived for the problem it posed. A longer result is accepted:
    /// it is readable, and the condition checked here is exactly the one that
    /// makes the reads legal.
    auto result_covers_problem(const qp_result<Scalar>& result) const -> bool
    {
        auto dims = qp_dimensions();
        auto const primal_size = static_cast<Eigen::Index>(dims.n_dec);
        auto const dual_size = static_cast<Eigen::Index>(dims.n_con);
        return result.x.size() >= primal_size
            && result.y.size() >= dual_size
            && result.x.head(primal_size).allFinite()
            && result.y.head(dual_size).allFinite()
            && std::isfinite(result.objective)
            && std::isfinite(result.solve_time)
            && std::isfinite(result.primal_residual)
            && std::isfinite(result.dual_residual)
            && result.iterations >= 0;
    }

    auto build_qp_structure(const Matrix<Scalar, NX, NX>& A_lin,
                            const Matrix<Scalar, NY, NX>& H_lin,
                            const cov_matrix_t& arrival_inverse)
        -> qp_problem<Scalar>
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        bool has_residual = m_residual_bound.has_value();
        std::array<Matrix<Scalar, NX, NX>, 1> A_arr{A_lin};
        std::array<Matrix<Scalar, NY, NX>, 1> H_arr{H_lin};

        return detail::build_mhe_qp_structure<Scalar, NX, NU, NY>(N, m_arrival_cost_weight, arrival_inverse, m_Q_inv, m_R_inv, A_arr, H_arr, has_box, m_soft_constraints && has_box, m_soft_penalty, has_residual);
    }

    auto build_qp_update(const Matrix<Scalar, NX, NX>& A_lin,
                         const Matrix<Scalar, NX, NU>& B_lin,
                         const Matrix<Scalar, NY, NX>& H_lin,
                         const cov_matrix_t& arrival_inverse)
        -> qp_update<Scalar>
    {
        bool has_box = m_x_min.has_value() || m_x_max.has_value();
        std::span<const input_vector_t> u_span{m_u_window.data(), N};
        std::span<const output_vector_t> z_span{m_z_window.data(), N + 1};

        return detail::build_mhe_qp_update<Scalar, NX, NU, NY>(N, m_arrival_cost_weight, arrival_inverse, m_Q_inv, m_R_inv, A_lin, B_lin, H_lin, m_prior_state_window[0], u_span, z_span, has_box, m_soft_constraints && has_box, m_soft_penalty, m_x_min, m_x_max, m_residual_bound, m_warm_z, m_warm_y);
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

    /// @brief Abandon the window solve and report the embedded filter's estimate
    /// instead. This is the estimator's whole failure channel: `update` returns
    /// nothing, so a caller learns what happened from `used_ekf_fallback` plus
    /// the reported status, and `reason` is what names the condition there.
    void fallback_to_ekf(solve_status reason = solve_status::error)
    {
        m_x_window[N] = m_ekf.state();
        m_innovation = m_ekf.innovation();
        m_diagnostics = mhe_diagnostics<Scalar>{.status = reason, .used_ekf_fallback = true};
    }

    // The slices below need no length check of their own: this runs only from
    // the extraction, which the solve attempt reaches only after confirming the
    // result covers the posed problem.
    void shift_warm_start(const Eigen::VectorX<Scalar>& sol_x, const Eigen::VectorX<Scalar>& sol_y)
    {
        for(int k = 0; k < Ni; ++k)
            m_warm_z.segment(k * nx, nx) = sol_x.segment((k + 1) * nx, nx);
        m_warm_z.segment(Ni * nx, nx) = m_ekf.state();

        shift_slack_warm_start(sol_x);
        m_warm_y = sol_y;
    }

    void shift_slack_warm_start(const Eigen::VectorX<Scalar>& /*sol_x*/)
    {
        auto dims = detail::compute_mhe_dims<NX, NY>(N, m_x_min.has_value() || m_x_max.has_value(), m_soft_constraints && (m_x_min.has_value() || m_x_max.has_value()), m_residual_bound.has_value());
        if(dims.n_slack <= 0)
            return;

        // Non-negative slacks warm-start at zero: the softened bound is inactive
        // unless the shifted estimate reenters the constraint region.
        m_warm_z.segment(dims.n_states, dims.n_slack).setZero();
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
    std::array<state_vector_t, N + 1> m_prior_state_window;
    std::array<cov_matrix_t, N + 1> m_prior_cov_window;

    std::size_t m_step_count{0};
    std::size_t m_update_count{0};
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

// This solver exists only to instantiate the observer-policy static assertions
// below. It has no setup failure mode, so its setup-error type carries no
// enumerators: the solver concept accepts one setup shape, a fallible one, and a
// backend with nothing to fail at writes a trivially succeeding fallible setup.
enum class mhe_sa_setup_error
{
};

struct mhe_sa_solver
{
    using scalar_type = double;

    auto setup(const qp_problem<double>&) -> ctrlpp::expected<void, mhe_sa_setup_error> { return {}; }

    auto solve(const qp_update<double>&) -> qp_result<double> { return {}; }
};

}

static_assert(ObserverPolicy<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);
static_assert(CovarianceObserver<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);

}

#endif
