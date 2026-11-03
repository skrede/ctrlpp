#ifndef HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H
#define HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H

// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960

#include "ctrlpp/estimation/kalman_types.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/estimation/detail/kalman_steps.h"
#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class kalman_filter : private detail::kalman_steps<Scalar, NX, NU, NY>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    using base = detail::kalman_steps<Scalar, NX, NU, NY>;

    struct validated_tag
    {
    };

public:
    using observer_tag = struct kalman_tag;
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;
    using typename base::cov_matrix_t;
    using typename base::meas_cov_matrix_t;
    using typename base::system_t;

    static ctrlpp::expected<kalman_filter, filter_error> create(system_t sys, kalman_config<Scalar, NX, NU, NY> config)
    {
        if(const auto valid = detail::validate_filter_configuration(config.Q, config.R, config.x0, config.P0); !valid)
            return ctrlpp::unexpected(valid.error());
        return kalman_filter{validated_tag{}, std::move(sys), std::move(config)};
    }

    // Not fallible: the plant already took the command. A prediction that poisons the estimate is latched
    // in health() by the next update.
    void predict(const input_vector_t& u)
    {
        this->m_P_post_prev = this->m_P;
        this->propagate_state(u);
        this->propagate_covariance();
    }

    // A refused step leaves the state, the covariance, the innovation and the NIS bitwise unchanged, so the
    // caller may retry with the next sample.
    ctrlpp::expected<void, kalman_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_x = this->m_x;
        auto const previous_P = this->m_P;
        auto const previous_innovation = this->m_innovation;
        auto const previous_nis = this->m_nis_value;
        this->compute_innovation(z);
        auto S = this->compute_innovation_covariance();
        auto K = this->compute_kalman_gain(S);

        this->apply_state_correction(K);
        this->update_covariance(K);
        this->compute_nis(S);
        if(!this->m_x.allFinite() || !this->m_P.allFinite() || !this->m_innovation.allFinite() || !std::isfinite(this->m_nis_value))
        {
            this->m_x = previous_x;
            this->m_P = previous_P;
            this->m_innovation = previous_innovation;
            this->m_nis_value = previous_nis;
            return ctrlpp::unexpected(kalman_update_error::non_finite_result);
        }
        return {};
    }

    const state_vector_t& state() const { return this->m_x; }

    const cov_matrix_t& covariance() const { return this->m_P; }

    const output_vector_t& innovation() const { return this->m_innovation; }

    Scalar nis() const { return this->m_nis_value; }

    kalman_health health() const { return this->m_health; }

    bool is_steady_state(Scalar tol = Scalar{1e-10}) const
    {
        Scalar p_norm = this->m_P.norm();
        if(p_norm < std::numeric_limits<Scalar>::epsilon())
            return true;
        return (this->m_P - this->m_P_post_prev).norm() / p_norm < tol;
    }

    // Refuses what create refuses, before assigning anything; a drift-asymmetric P0 is stored symmetrized.
    ctrlpp::expected<void, filter_error> reset_covariance(cov_matrix_t P0)
    {
        if(const auto valid = detail::validate_initial_covariance(P0); !valid)
            return ctrlpp::unexpected(valid.error());
        this->m_P = P0;
        this->m_P_post_prev = P0;
        return {};
    }

    void set_model(system_t sys) { this->m_sys = std::move(sys); }

    // Assigns both matrices or neither; a drift-asymmetric matrix is stored symmetrized.
    ctrlpp::expected<void, filter_error> set_noise(cov_matrix_t Q, meas_cov_matrix_t R)
    {
        if(const auto valid = detail::validate_noise(Q, R); !valid)
            return ctrlpp::unexpected(valid.error());
        this->m_Q = std::move(Q);
        this->m_R = std::move(R);
        return {};
    }

private:
    kalman_filter(validated_tag, system_t sys, kalman_config<Scalar, NX, NU, NY> config)
        : base{std::move(sys), std::move(config)}
    {
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, kalman_update_error> check_step(const output_vector_t& z) const
    {
        if(!this->m_x.allFinite())
            return ctrlpp::unexpected(kalman_update_error::non_finite_state);
        if(!this->m_P.allFinite())
            return ctrlpp::unexpected(kalman_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(kalman_update_error::non_finite_measurement);
        return {};
    }

    kalman_update_error latch_health(kalman_update_error fault)
    {
        if(fault != kalman_update_error::non_finite_measurement)
            this->m_health = kalman_health::non_finite_estimate;
        return fault;
    }
};

static_assert(ObserverPolicy<kalman_filter<double, 2, 1, 1>>);
static_assert(CovarianceObserver<kalman_filter<double, 2, 1, 1>>);

}

#endif
