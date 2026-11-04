#ifndef HPP_GUARD_CTRLPP_MHE_H
#define HPP_GUARD_CTRLPP_MHE_H

// C. V. Rao, J. B. Rawlings, J. H. Lee, "Constrained linear state estimation -- a moving horizon approach", Automatica 37(10), 2001.
// @cite rao2003 -- Rao, Rawlings & Mayne, "Constrained State Estimation for Nonlinear Discrete-Time Systems", 2003
// @cite kuhl2011 -- Kuhl et al., "A Real-Time Algorithm for Moving Horizon State and Parameter Estimation", 2011

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"

#include "ctrlpp/mhe/detail/mhe_solve.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include <span>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, typename Solver, typename Dynamics, typename Measurement>
    requires qp_solver<Solver> && dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class mhe : private detail::mhe_solve<Scalar, NX, NU, NY, N, Solver, Dynamics, Measurement>
{
    static_assert(N > 0,
            "Window length N must be positive: it sizes the fixed estimation window arrays, which the update rotates, reads the trailing element of, and indexes at their midpoint, "
            "none of which is defined for an empty window");

    using base = detail::mhe_solve<Scalar, NX, NU, NY, N, Solver, Dynamics, Measurement>;

    // Selects the non-validating constructor reserved for create, so the factory's validation cannot be bypassed.
    struct validated_tag
    {
    };

public:
    using observer_tag = struct mhe_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;

    static ctrlpp::expected<mhe, moving_horizon_construction_error> create(Dynamics dynamics, Measurement measurement, const mhe_config<Scalar, NX, NU, NY, N>& config)
    {
        auto filter = ekf<Scalar, NX, NU, NY, Dynamics, Measurement>::create(
                dynamics, measurement, ekf_config<Scalar, NX, NU, NY>{.Q = config.Q, .R = config.R, .x0 = config.x0, .P0 = config.P0, .numerical_eps = config.numerical_eps});
        if(!filter)
            return ctrlpp::unexpected(moving_horizon_construction_error{filter.error()});

        auto q_inv = detail::finite_full_piv_inverse(config.Q, moving_horizon_configuration_error::non_invertible_process_noise);
        if(!q_inv)
            return ctrlpp::unexpected(moving_horizon_construction_error{q_inv.error()});
        auto r_inv = detail::finite_full_piv_inverse(config.R, moving_horizon_configuration_error::non_invertible_measurement_noise);
        if(!r_inv)
            return ctrlpp::unexpected(moving_horizon_construction_error{r_inv.error()});
        auto p0_inv = detail::finite_full_piv_inverse(config.P0, moving_horizon_configuration_error::non_invertible_initial_covariance);
        if(!p0_inv)
            return ctrlpp::unexpected(moving_horizon_construction_error{p0_inv.error()});
        auto options = detail::validate_moving_horizon_options(config);
        if(!options)
            return ctrlpp::unexpected(moving_horizon_construction_error{options.error()});

        return mhe{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(*filter), config, std::move(*q_inv), std::move(*r_inv)};
    }

    void predict(const input_vector_t& u)
    {
        this->m_ekf.predict(u);
        std::rotate(this->m_u_window.begin(), this->m_u_window.begin() + 1, this->m_u_window.end());
        this->m_u_window.back() = u;
        ++this->m_step_count;
    }

    ctrlpp::expected<void, ekf_update_error> update(const output_vector_t& z)
    {
        // The prior is captured before the correction: shifted into the window, the oldest entry is then the estimate
        // at the window head from data strictly before the window, the arrival prior that does not double-count.
        const state_vector_t prior_state = this->m_ekf.state();
        const cov_matrix_t prior_cov = this->m_ekf.covariance();

        if(const auto stepped = this->m_ekf.update(z); !stepped)
        {
            this->invalidate_window_after_refusal();
            return ctrlpp::unexpected(stepped.error());
        }

        this->shift_windows(prior_state, prior_cov, z);
        if(this->m_update_count <= N)
        {
            this->use_filter_estimate();
            return {};
        }

        this->solve_mhe(z);
        return {};
    }

    const state_vector_t& state() const { return this->m_x_window[N]; }
    const cov_matrix_t& covariance() const { return this->m_ekf.covariance(); }
    const output_vector_t& innovation() const { return this->m_innovation; }
    std::span<const state_vector_t> trajectory() const { return {this->m_x_window.data(), N + 1}; }
    const state_vector_t& arrival_state() const { return this->m_x_window[0]; }
    const cov_matrix_t& arrival_covariance() const { return this->m_prior_cov_window[0]; }
    bool is_initialized() const { return this->m_update_count > N; }
    const mhe_diagnostics<Scalar>& diagnostics() const { return this->m_diagnostics; }

private:
    mhe(validated_tag, Dynamics dynamics, Measurement measurement, ekf<Scalar, NX, NU, NY, Dynamics, Measurement> filter, const mhe_config<Scalar, NX, NU, NY, N>& config,
            cov_matrix_t q_inv, Matrix<Scalar, NY, NY> r_inv)
        : base{std::move(dynamics), std::move(measurement), std::move(filter), config, std::move(q_inv), std::move(r_inv)}
    {
    }
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

// Exists only to instantiate the observer-policy assertions below. The solver concept accepts only a fallible setup,
// so a backend with nothing to fail at declares an error type without enumerators.
enum class mhe_sa_setup_error
{
};

struct mhe_sa_solver
{
    using scalar_type = double;

    ctrlpp::expected<void, mhe_sa_setup_error> setup(const qp_problem<double>&) { return {}; }

    qp_result<double> solve(const qp_update<double>&) { return {}; }
};

}

static_assert(ObserverPolicy<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);
static_assert(CovarianceObserver<mhe<double, 2, 1, 1, 5, detail::mhe_sa_solver, detail::mhe_sa_dynamics, detail::mhe_sa_measurement>>);

}

#endif
