#ifndef HPP_GUARD_CTRLPP_ESTIMATION_UKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_UKF_H

// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001

#include "ctrlpp/estimation/ukf_types.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include "ctrlpp/estimation/detail/ukf_steps.h"
#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy = merwe_sigma_points<Scalar, NX>>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && sigma_point_strategy<Strategy, Scalar, NX>
class ukf : private detail::ukf_steps<Scalar, NX, NU, NY, Dynamics, Measurement, Strategy>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    using base = detail::ukf_steps<Scalar, NX, NU, NY, Dynamics, Measurement, Strategy>;

    struct validated_tag
    {
    };

public:
    using observer_tag = struct ukf_tag;
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;
    using typename base::cov_matrix_t;
    using typename base::meas_cov_matrix_t;

    static ctrlpp::expected<ukf, filter_error> create(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, Strategy strategy = Strategy{})
    {
        if(const auto valid = detail::validate_filter_configuration(config.Q, config.R, config.x0, config.P0); !valid)
            return ctrlpp::unexpected(valid.error());
        return ukf{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config), std::move(strategy)};
    }

    static ctrlpp::expected<ukf, filter_error> create(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, typename Strategy::options_t strategy_options)
    {
        auto strategy = Strategy::try_create(std::move(strategy_options));
        if(!strategy)
            return ctrlpp::unexpected(strategy.error());
        return create(std::move(dynamics), std::move(measurement), std::move(config), std::move(*strategy));
    }

    // Not fallible: the plant already took the command. A prediction that poisons the estimate or is not
    // positive definite is latched in health() instead.
    void predict(const input_vector_t& u)
    {
        auto sigma = this->m_strategy.generate(this->m_x, this->m_P);
        this->note_health(sigma);
        auto propagated = this->propagate_sigma_points(sigma.points, u);
        auto x_pred = this->compute_predicted_mean(sigma.Wm, propagated);
        this->m_P = this->compute_predicted_covariance(sigma.Wc, propagated, x_pred);
        this->m_x = x_pred;
        this->note_prediction();
    }

    // A refused step leaves every member bitwise unchanged, the health included, so the caller may retry
    // with the next sample.
    ctrlpp::expected<void, ukf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_x = this->m_x;
        auto const previous_P = this->m_P;
        auto const previous_innovation = this->m_innovation;
        auto const previous_nis = this->m_nis;
        auto const previous_health = this->m_health;
        auto sigma = this->m_strategy.generate(this->m_x, this->m_P);
        this->note_health(sigma);
        auto z_sigma = this->compute_measurement_sigma_points(sigma.points);
        auto z_pred = this->compute_predicted_measurement(sigma.Wm, z_sigma);
        auto [S, Pxz] = this->compute_innovation_and_cross_covariance(sigma, z_sigma, z_pred);
        auto K = this->compute_kalman_gain(Pxz, S);

        this->m_innovation = (z - z_pred).eval();
        this->m_nis = (this->m_innovation.transpose() * S.colPivHouseholderQr().solve(this->m_innovation))(0, 0);
        this->apply_correction_and_update_covariance(K, S, this->reduction_base(sigma));
        if(!this->m_x.allFinite() || !this->m_P.allFinite() || !this->m_innovation.allFinite() || !std::isfinite(this->m_nis))
            return roll_back(previous_x, previous_P, previous_innovation, previous_nis, previous_health);
        this->note_definiteness();
        return {};
    }

    const state_vector_t& state() const { return this->m_x; }

    const cov_matrix_t& covariance() const { return this->m_P; }

    const output_vector_t& innovation() const { return this->m_innovation; }

    Scalar nis() const { return this->m_nis; }

    ukf_health health() const { return this->m_health; }

private:
    ukf(validated_tag, Dynamics&& dynamics, Measurement&& measurement, ukf_config<Scalar, NX, NU, NY>&& config, Strategy&& strategy)
        : base{std::move(dynamics), std::move(measurement), std::move(config), std::move(strategy)}
    {
    }

    ctrlpp::expected<void, ukf_update_error> roll_back(const state_vector_t& x, const cov_matrix_t& P, const output_vector_t& innovation, Scalar nis, ukf_health health)
    {
        this->m_x = x;
        this->m_P = P;
        this->m_innovation = innovation;
        this->m_nis = nis;
        this->m_health = health;
        return ctrlpp::unexpected(ukf_update_error::non_finite_result);
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, ukf_update_error> check_step(const output_vector_t& z) const
    {
        if(!this->m_x.allFinite())
            return ctrlpp::unexpected(ukf_update_error::non_finite_state);
        if(!this->m_P.allFinite())
            return ctrlpp::unexpected(ukf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(ukf_update_error::non_finite_measurement);
        return {};
    }

    ukf_update_error latch_health(ukf_update_error fault)
    {
        if(fault != ukf_update_error::non_finite_measurement)
            this->m_health = ukf_health::non_finite_estimate;
        return fault;
    }
};

template <typename Dynamics, typename Measurement, ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
ukf(Dynamics, Measurement, ukf_config<Scalar, NX, NU, NY>) -> ukf<Scalar, NX, NU, NY, Dynamics, Measurement, merwe_sigma_points<Scalar, NX>>;

namespace detail
{

struct ukf_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>&, const Vector<double, 1>&) const { return Vector<double, 2>::Zero(); }
};

struct ukf_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>&) const { return Vector<double, 1>::Zero(); }
};

}

static_assert(ObserverPolicy<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);

}

#endif
