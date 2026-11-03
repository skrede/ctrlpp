#ifndef HPP_GUARD_CTRLPP_ESTIMATION_EKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_EKF_H

// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13

#include "ctrlpp/estimation/ekf_types.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/estimation/detail/ekf_steps.h"
#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <cmath>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class ekf : private detail::ekf_steps<Scalar, NX, NU, NY, Dynamics, Measurement>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    using base = detail::ekf_steps<Scalar, NX, NU, NY, Dynamics, Measurement>;

    struct validated_tag
    {
    };

public:
    using observer_tag = struct ekf_tag;
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;
    using typename base::cov_matrix_t;
    using typename base::meas_cov_matrix_t;

    // numerical_eps is not validated: it is read only when a model is not analytically differentiable.
    static ctrlpp::expected<ekf, filter_error> create(Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config)
    {
        if(const auto valid = detail::validate_filter_configuration(config.Q, config.R, config.x0, config.P0); !valid)
            return ctrlpp::unexpected(valid.error());
        return ekf{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config)};
    }

    // Not fallible: the plant already took the command. A prediction that poisons the estimate is latched
    // in health() by the next update.
    void predict(const input_vector_t& u)
    {
        auto x_prev = this->m_x;
        this->propagate_state(x_prev, u);
        this->propagate_covariance(x_prev, u);
    }

    // A refused step leaves the state, the covariance, the innovation and the NIS bitwise unchanged, so the
    // caller may retry with the next sample.
    ctrlpp::expected<void, ekf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_x = this->m_x;
        auto const previous_P = this->m_P;
        auto const previous_innovation = this->m_innovation;
        auto const previous_nis = this->m_nis;
        auto z_pred = this->m_measurement(this->m_x);
        this->m_innovation = (z - z_pred).eval();

        auto H = this->compute_measurement_jacobian();
        auto S = this->compute_innovation_covariance(H);
        auto K = this->compute_kalman_gain(H, S);

        this->apply_state_correction(K, z, z_pred);
        this->update_covariance(K, H);

        this->m_nis = (this->m_innovation.transpose() * S.colPivHouseholderQr().solve(this->m_innovation))(0, 0);
        if(!this->m_x.allFinite() || !this->m_P.allFinite() || !this->m_innovation.allFinite() || !std::isfinite(this->m_nis))
        {
            this->m_x = previous_x;
            this->m_P = previous_P;
            this->m_innovation = previous_innovation;
            this->m_nis = previous_nis;
            return ctrlpp::unexpected(ekf_update_error::non_finite_result);
        }
        return {};
    }

    const state_vector_t& state() const { return this->m_x; }

    const cov_matrix_t& covariance() const { return this->m_P; }

    const output_vector_t& innovation() const { return this->m_innovation; }

    Scalar nis() const { return this->m_nis; }

    ekf_health health() const { return this->m_health; }

private:
    ekf(validated_tag, Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config)
        : base{std::move(dynamics), std::move(measurement), std::move(config)}
    {
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, ekf_update_error> check_step(const output_vector_t& z) const
    {
        if(!this->m_x.allFinite())
            return ctrlpp::unexpected(ekf_update_error::non_finite_state);
        if(!this->m_P.allFinite())
            return ctrlpp::unexpected(ekf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(ekf_update_error::non_finite_measurement);
        return {};
    }

    ekf_update_error latch_health(ekf_update_error fault)
    {
        if(fault != ekf_update_error::non_finite_measurement)
            this->m_health = ekf_health::non_finite_estimate;
        return fault;
    }
};

namespace detail
{

struct ekf_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>&, const Vector<double, 1>&) const { return Vector<double, 2>::Zero(); }
};

struct ekf_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>&) const { return Vector<double, 1>::Zero(); }
};

}

static_assert(ObserverPolicy<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);
static_assert(CovarianceObserver<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);

}

#endif
