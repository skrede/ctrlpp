#ifndef HPP_GUARD_CTRLPP_NMHE_H
#define HPP_GUARD_CTRLPP_NMHE_H

// @cite diehl2009 -- Diehl, Ferreau & Haverbeke, "Efficient Numerical Methods for Nonlinear MPC and Moving Horizon Estimation", 2009
// V. M. Zavala, L. T. Biegler, "The advanced-step NMPC controller: Optimality, stability and robustness", Automatica 45(1), 2009.
// @cite rao2003 -- Rao, Rawlings & Mayne, "Constrained State Estimation for Nonlinear Discrete-Time Systems", 2003

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/mhe/mhe_config.h"
#include "ctrlpp/mhe/mhe_diagnostics.h"
#include "ctrlpp/mhe/mhe_nlp_formulation.h"

#include "ctrlpp/mhe/detail/nmhe_solve.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include <span>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, nlp_solver Solver, typename Dynamics, typename Measurement, std::size_t NC = 0>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class nmhe : private detail::nmhe_solve<Scalar, NX, NU, NY, N, NC, Solver, Dynamics, Measurement>
{
    static_assert(N > 0,
            "Window length N must be positive: it sizes the fixed estimation window arrays, which the update rotates and reads the trailing element of, neither of which is defined for "
            "an empty window");

    using base = detail::nmhe_solve<Scalar, NX, NU, NY, N, NC, Solver, Dynamics, Measurement>;

    // Selects the non-validating constructor reserved for create, so the factory's validation cannot be bypassed.
    struct validated_tag
    {
    };

public:
    using observer_tag = struct nmhe_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;

    static ctrlpp::expected<nmhe, moving_horizon_construction_error> create(Dynamics dynamics, Measurement measurement, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config)
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
        auto path = validate_path_constraint(config);
        if(!path)
            return ctrlpp::unexpected(moving_horizon_construction_error{path.error()});

        return nmhe{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(*filter), config, std::move(*q_inv), std::move(*r_inv), std::move(*p0_inv)};
    }

    void predict(const input_vector_t& u)
    {
        this->m_ekf.predict(u);
        std::rotate(this->m_u_window.begin(), this->m_u_window.begin() + 1, this->m_u_window.end());
        this->m_u_window.back() = u;
        this->m_state->u_window = this->m_u_window;
        ++this->m_step_count;
    }

    // The prior is captured before the correction: shifted into the window, the oldest entry is then the estimate at the
    // window head from data strictly before the window, the arrival prior that does not double-count.
    ctrlpp::expected<void, ekf_update_error> update(const output_vector_t& z)
    {
        const state_vector_t prior_state = this->m_ekf.state();
        const cov_matrix_t prior_cov = this->m_ekf.covariance();

        if(const auto stepped = this->m_ekf.update(z); !stepped)
        {
            this->invalidate_window_after_refusal();
            return ctrlpp::unexpected(stepped.error());
        }

        this->shift_priors(prior_state, prior_cov);
        this->shift_measurements(z);
        if(this->m_update_count <= N)
        {
            if(this->m_step_count < N)
                this->m_x_window[this->m_step_count] = this->m_ekf.state();
            this->m_x_window[N] = this->m_ekf.state();
            this->m_innovation = this->m_ekf.innovation();
            this->m_diagnostics = detail::ekf_fallback_diagnostics<Scalar>(solve_status::optimal);
            return {};
        }

        this->solve_nmhe();
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
    nmhe(validated_tag, Dynamics dynamics, Measurement measurement, ekf<Scalar, NX, NU, NY, Dynamics, Measurement> filter, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config,
            cov_matrix_t q_inv, Matrix<Scalar, NY, NY> r_inv, cov_matrix_t p0_inv)
        : base{std::move(dynamics), std::move(measurement), std::move(filter), config, std::move(q_inv), std::move(r_inv), std::move(p0_inv)}
    {
    }

    static ctrlpp::expected<void, moving_horizon_configuration_error> validate_path_constraint(const nmhe_config<Scalar, NX, NU, NY, N, NC>& config)
    {
        if constexpr(NC > 0)
        {
            if(config.path_constraint_cb)
            {
                if(!config.path_penalty.allFinite() || !(config.path_penalty.array() > Scalar{0}).all())
                    return ctrlpp::unexpected(moving_horizon_configuration_error::invalid_path_penalty);
                if(!(*config.path_constraint_cb)(config.x0).allFinite())
                    return ctrlpp::unexpected(moving_horizon_configuration_error::non_finite_path_constraint);
            }
        }
        return {};
    }
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

// Exists only to instantiate the observer-policy assertions below. The solver concept accepts only a fallible setup,
// so a backend with nothing to fail at declares an error type without enumerators.
enum class nmhe_sa_setup_error
{
};

struct nmhe_sa_solver
{
    using scalar_type = double;

    ctrlpp::expected<void, nmhe_sa_setup_error> setup(const nlp_problem<double>&) { return {}; }

    nlp_result<double> solve(const nlp_update<double>&) { return {}; }
};

}

static_assert(ObserverPolicy<nmhe<double, 2, 1, 1, 5, detail::nmhe_sa_solver, detail::nmhe_sa_dynamics, detail::nmhe_sa_measurement>>);
static_assert(CovarianceObserver<nmhe<double, 2, 1, 1, 5, detail::nmhe_sa_solver, detail::nmhe_sa_dynamics, detail::nmhe_sa_measurement>>);

}

#endif
