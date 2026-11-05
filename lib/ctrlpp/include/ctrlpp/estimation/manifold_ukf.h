#ifndef HPP_GUARD_CTRLPP_ESTIMATION_MANIFOLD_UKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_MANIFOLD_UKF_H

// @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/manifold_ukf_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/estimation/sigma_points/so3_sigma_points.h"

#include "ctrlpp/estimation/detail/manifold_ukf_steps.h"
#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <Eigen/Cholesky>
#include <Eigen/Geometry>

#include <cmath>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy = so3_merwe_sigma_points<Scalar>>
    requires manifold_ukf_dynamics_model<Dynamics, Scalar> && manifold_ukf_measurement_model<Measurement, Scalar, NY> && manifold_sigma_point_strategy<Strategy, Scalar>
class manifold_ukf : private detail::manifold_ukf_steps<Scalar, NY, Dynamics, Measurement, Strategy>
{
    static_assert(NY > 0, "Output dimension NY must be positive");

    using base = detail::manifold_ukf_steps<Scalar, NY, Dynamics, Measurement, Strategy>;

    struct validated_tag
    {
    };

public:
    using observer_tag = struct manifold_ukf_tag;
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;
    using typename base::cov_matrix_t;
    using typename base::meas_cov_t;

    static ctrlpp::expected<manifold_ukf, filter_error> create(Dynamics dynamics, Measurement measurement, manifold_ukf_config<Scalar, NY> config, Strategy strategy = Strategy{})
    {
        if(const auto attitude = detail::checked_unit_quaternion(config.q0); !attitude)
            return ctrlpp::unexpected(attitude.error());
        if(const auto noise = detail::validate_noise(config.Q, config.R); !noise)
            return ctrlpp::unexpected(noise.error());
        if(const auto prior = detail::validate_initial_covariance(config.P0); !prior)
            return ctrlpp::unexpected(prior.error());
        return manifold_ukf{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config), std::move(strategy)};
    }

    // Not fallible: the body already performed the rotation. A prediction that poisons the attitude or is
    // not positive definite is latched in health() instead.
    void predict(const input_vector_t& omega)
    {
        auto sigma = this->m_strategy.generate(this->m_q, this->m_P);
        note_health(sigma);
        auto q_prop = this->propagate_manifold_sigma_points(sigma.points, omega);
        auto mean_result = this->compute_manifold_predicted_mean(q_prop, sigma.Wm);
        auto q_new = mean_result.mean;
        this->m_P = this->compute_manifold_predicted_covariance(q_prop, q_new, sigma.Wc);
        this->m_q = q_new;
        this->update_state_cache();
        note_prediction(mean_result.converged);
    }

    // A refused step leaves every member bitwise unchanged, so the caller may retry with the next sample.
    ctrlpp::expected<void, manifold_ukf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = this->m_q;
        auto const previous_P = this->m_P;
        auto const previous_innovation = this->m_innovation;
        auto const previous_state = this->m_state_cache;
        auto sigma = this->m_strategy.generate(this->m_q, this->m_P);
        auto z_sigma = this->compute_measurement_sigma_points(sigma.points);
        auto z_pred = this->compute_predicted_measurement(sigma.Wm, z_sigma);
        auto [S, Pxz] = this->compute_innovation_and_cross_covariance(sigma, z_sigma, z_pred);
        auto K = this->compute_kalman_gain(Pxz, S);

        this->m_innovation = (z - z_pred).eval();
        this->apply_manifold_correction(K, S, this->reduction_base(sigma));
        this->update_state_cache();
        if(!this->m_q.coeffs().allFinite() || !this->m_P.allFinite() || !this->m_innovation.allFinite() || !this->m_state_cache.allFinite())
            return roll_back(previous_q, previous_P, previous_innovation, previous_state);
        note_health(sigma);
        note_definiteness();
        return {};
    }

    const state_vector_t& state() const { return this->m_state_cache; }

    const cov_matrix_t& covariance() const { return this->m_P; }

    const output_vector_t& innovation() const { return this->m_innovation; }

    const Eigen::Quaternion<Scalar>& attitude() const { return this->m_q; }

    manifold_ukf_health health() const { return this->m_health; }

private:
    manifold_ukf(validated_tag, Dynamics dynamics, Measurement measurement, manifold_ukf_config<Scalar, NY> config, Strategy strategy)
        : base{std::move(dynamics), std::move(measurement), std::move(config), std::move(strategy)}
    {
    }

    ctrlpp::expected<void, manifold_ukf_update_error> roll_back(
            const Eigen::Quaternion<Scalar>& q, const cov_matrix_t& P, const output_vector_t& innovation, const state_vector_t& state)
    {
        this->m_q = q;
        this->m_P = P;
        this->m_innovation = innovation;
        this->m_state_cache = state;
        return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_result);
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, manifold_ukf_update_error> check_step(const output_vector_t& z) const
    {
        if(!this->m_q.coeffs().allFinite())
            return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_state);
        if(!this->m_P.allFinite())
            return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_measurement);
        return {};
    }

    manifold_ukf_update_error latch_health(manifold_ukf_update_error fault)
    {
        if(fault != manifold_ukf_update_error::non_finite_measurement)
            this->m_health = manifold_ukf_health::non_finite_estimate;
        return fault;
    }

    template <typename Sigma>
    void note_health(const Sigma& sigma)
    {
        if(sigma.spd_repaired && this->m_health == manifold_ukf_health::ok)
            this->m_health = manifold_ukf_health::covariance_repaired;
    }

    void note_prediction(bool mean_converged)
    {
        if(!this->m_q.coeffs().allFinite() || !this->m_P.allFinite())
            this->m_health = manifold_ukf_health::non_finite_estimate;
        else if(!mean_converged && this->m_health < manifold_ukf_health::mean_not_converged)
            this->m_health = manifold_ukf_health::mean_not_converged;
        else
            note_definiteness();
    }

    // Judged by the unpivoted Cholesky factorization the next sigma-point generation attempts first on this
    // same matrix, so the status never disagrees with the repair that generation makes.
    void note_definiteness()
    {
        if(this->m_health == manifold_ukf_health::ok && Eigen::LLT<cov_matrix_t>{this->m_P}.info() != Eigen::Success)
            this->m_health = manifold_ukf_health::covariance_repaired;
    }
};

namespace detail
{

struct mukf_sa_dynamics
{
    Eigen::Quaternion<double> operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>&) const { return q; }
};

struct mukf_sa_measurement
{
    Vector<double, 3> operator()(const Eigen::Quaternion<double>&) const { return Vector<double, 3>::Zero(); }
};

}

static_assert(ObserverPolicy<manifold_ukf<double, 3, detail::mukf_sa_dynamics, detail::mukf_sa_measurement>>);
static_assert(CovarianceObserver<manifold_ukf<double, 3, detail::mukf_sa_dynamics, detail::mukf_sa_measurement>>);

}

#endif
