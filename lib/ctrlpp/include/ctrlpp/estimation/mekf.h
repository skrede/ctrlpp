#ifndef HPP_GUARD_CTRLPP_ESTIMATION_MEKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_MEKF_H

// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
// @cite crassidis2003 -- Crassidis & Markley, "Unscented Filtering for Spacecraft Attitude Estimation", J. Guidance Control Dyn 26(4), 2003

#include "ctrlpp/estimation/mekf_types.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/estimation/detail/mekf_steps.h"
#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <Eigen/Geometry>

#include <cmath>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NB, std::size_t NY, typename Measurement>
    requires mekf_measurement_model<Measurement, Scalar, NB, NY>
class mekf : private detail::mekf_steps<Scalar, NB, NY, Measurement>
{
    static_assert(NB >= 3, "Bias dimension NB must be at least 3: the propagation subtracts the leading three bias elements from the gyro rate");
    static_assert(NY > 0, "Output dimension NY must be positive");

    using base = detail::mekf_steps<Scalar, NB, NY, Measurement>;

    struct validated_tag
    {
    };

public:
    using observer_tag = struct mekf_tag;
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;
    using typename base::cov_matrix_t;
    using typename base::meas_cov_t;

    // q0 is checked before the normalization that seeds the state, since a zero or non-finite norm would
    // normalize to NaN.
    static ctrlpp::expected<mekf, filter_error> create(Measurement measurement, mekf_config<Scalar, NB, NY> config)
    {
        if(const auto attitude = detail::checked_unit_quaternion(config.q0); !attitude)
            return ctrlpp::unexpected(attitude.error());
        if(const auto noise = detail::validate_noise(config.Q, config.R); !noise)
            return ctrlpp::unexpected(noise.error());
        if(!config.b0.allFinite())
            return ctrlpp::unexpected(filter_error::non_finite_initial_bias);
        if(const auto prior = detail::validate_initial_covariance(config.P0); !prior)
            return ctrlpp::unexpected(prior.error());
        return mekf{validated_tag{}, std::move(measurement), std::move(config)};
    }

    // Not fallible: the body already rotated. A prediction that poisons the estimate is latched in
    // health() by the next update.
    void predict(const input_vector_t& omega) { this->predict_impl(omega, this->m_dt); }

    void predict(const input_vector_t& omega, Scalar dt) { this->predict_impl(omega, dt); }

    // A refused step leaves the attitude, the bias, the covariance and the innovation bitwise unchanged, so
    // the caller may retry with the next sample.
    ctrlpp::expected<void, mekf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = this->m_q;
        auto const previous_b = this->m_b;
        auto const previous_P = this->m_P;
        auto const previous_innovation = this->m_innovation;
        auto const previous_state = this->m_state_cache;
        auto z_pred = this->m_measurement(this->m_q, this->m_b);
        this->m_innovation = (z - z_pred).eval();

        auto H = this->compute_measurement_jacobian();
        auto S = this->compute_innovation_covariance(H);
        auto K = this->compute_kalman_gain(H, S);

        auto delta_xi = this->apply_multiplicative_correction(K);
        this->update_covariance(K, H, delta_xi);
        this->update_state_cache();
        if(!this->m_q.coeffs().allFinite() || !this->m_b.allFinite() || !this->m_P.allFinite() || !this->m_innovation.allFinite() || !this->m_state_cache.allFinite())
        {
            this->m_q = previous_q;
            this->m_b = previous_b;
            this->m_P = previous_P;
            this->m_innovation = previous_innovation;
            this->m_state_cache = previous_state;
            return ctrlpp::unexpected(mekf_update_error::non_finite_result);
        }
        return {};
    }

    const state_vector_t& state() const { return this->m_state_cache; }

    const cov_matrix_t& covariance() const { return this->m_P; }

    const output_vector_t& innovation() const { return this->m_innovation; }

    Eigen::Quaternion<Scalar> attitude() const { return this->m_q; }

    const Vector<Scalar, NB>& bias() const { return this->m_b; }

    mekf_health health() const { return this->m_health; }

private:
    mekf(validated_tag, Measurement measurement, mekf_config<Scalar, NB, NY> config)
        : base{std::move(measurement), std::move(config)}
    {
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, mekf_update_error> check_step(const output_vector_t& z) const
    {
        if(!this->m_q.coeffs().allFinite() || !this->m_b.allFinite())
            return ctrlpp::unexpected(mekf_update_error::non_finite_state);
        if(!this->m_P.allFinite())
            return ctrlpp::unexpected(mekf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(mekf_update_error::non_finite_measurement);
        return {};
    }

    mekf_update_error latch_health(mekf_update_error fault)
    {
        if(fault != mekf_update_error::non_finite_measurement)
            this->m_health = mekf_health::non_finite_estimate;
        return fault;
    }
};

namespace detail
{

struct mekf_sa_measurement
{
    Vector<double, 3> operator()(const Eigen::Quaternion<double>&, const Vector<double, 3>&) const { return Vector<double, 3>::Zero(); }
};

}

static_assert(ObserverPolicy<mekf<double, 3, 3, detail::mekf_sa_measurement>>);
static_assert(CovarianceObserver<mekf<double, 3, 3, detail::mekf_sa_measurement>>);

}

#endif
