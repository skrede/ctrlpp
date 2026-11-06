#ifndef HPP_GUARD_CTRLPP_ESTIMATION_COMPLEMENTARY_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_COMPLEMENTARY_FILTER_H

// @cite mahony2008 -- Mahony et al., "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/estimation/detail/filter_configuration.h"
#include "ctrlpp/estimation/detail/cf_attitude_correction.h"

#include <Eigen/Geometry>

#include <cmath>
#include <utility>

namespace ctrlpp
{

enum class cf_update_error
{
    non_finite_state,
    non_finite_measurement,
    non_finite_timestep,
    non_finite_result,
};

enum class cf_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar>
struct cf_config
{
    Scalar k_p{Scalar{2}};
    Scalar k_i{static_cast<Scalar>(0.005)};
    Scalar dt{static_cast<Scalar>(0.01)};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
};

template <ctrlpp_floating_scalar Scalar>
class complementary_filter
{
    struct validated_tag
    {
    };

public:
    using observer_tag = struct complementary_filter_tag;
    using state_vector_t = Vector<Scalar, 7>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, 3>;

    static ctrlpp::expected<complementary_filter, filter_error> create(cf_config<Scalar> config)
    {
        if(const auto attitude = detail::checked_unit_quaternion(config.q0); !attitude)
            return ctrlpp::unexpected(attitude.error());
        return complementary_filter{validated_tag{}, std::move(config)};
    }

    // A reading with no direction drops only its own correction term; the rate is integrated either way.
    // @cite mahony2008 -- Mahony, Hamel & Pflimlin, "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008, Sec. V-A
    ctrlpp::expected<void, cf_update_error> update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt)
    {
        if(const auto step = check_step(gyro, accel, dt); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = m_q;
        auto const previous_bias = m_bias;
        auto const previous_state = m_state_cache;
        Vector<Scalar, 3> direction;
        Vector<Scalar, 3> correction = Vector<Scalar, 3>::Zero();
        if(detail::unit_direction(accel, direction))
            correction = detail::compute_gravity_correction(m_q, direction);

        integrate_gyro(gyro, correction, dt);
        if(!m_q.coeffs().allFinite() || !m_bias.allFinite() || !m_state_cache.allFinite())
        {
            m_q = previous_q;
            m_bias = previous_bias;
            m_state_cache = previous_state;
            return ctrlpp::unexpected(cf_update_error::non_finite_result);
        }
        return {};
    }

    // @cite mahony2008 -- Mahony, Hamel & Pflimlin, "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008, Sec. V, Eq. 32
    ctrlpp::expected<void, cf_update_error> update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, const Vector<Scalar, 3>& mag, Scalar dt)
    {
        if(const auto step = check_step(gyro, accel, dt, &mag); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = m_q;
        auto const previous_bias = m_bias;
        auto const previous_state = m_state_cache;
        Vector<Scalar, 3> direction;
        Vector<Scalar, 3> correction = Vector<Scalar, 3>::Zero();
        if(detail::unit_direction(accel, direction))
            correction += detail::compute_gravity_correction(m_q, direction);
        if(detail::unit_direction(mag, direction))
            correction += detail::compute_magnetic_correction(m_q, direction);

        integrate_gyro(gyro, correction, dt);
        if(!m_q.coeffs().allFinite() || !m_bias.allFinite() || !m_state_cache.allFinite())
        {
            m_q = previous_q;
            m_bias = previous_bias;
            m_state_cache = previous_state;
            return ctrlpp::unexpected(cf_update_error::non_finite_result);
        }
        return {};
    }

    void predict(const input_vector_t& u) { m_gyro_buf = u; }

    ctrlpp::expected<void, cf_update_error> update(const output_vector_t& z) { return update(m_gyro_buf, z, m_dt); }

    const state_vector_t& state() const { return m_state_cache; }
    Eigen::Quaternion<Scalar> attitude() const { return m_q; }
    const Vector<Scalar, 3>& bias() const { return m_bias; }

    cf_health health() const { return m_health; }

private:
    Eigen::Quaternion<Scalar> m_q;
    Vector<Scalar, 3> m_bias;
    Scalar m_k_p;
    Scalar m_k_i;
    Scalar m_dt;
    Vector<Scalar, 3> m_gyro_buf;
    state_vector_t m_state_cache;
    cf_health m_health;

    complementary_filter(validated_tag, cf_config<Scalar> config)
        : m_q{config.q0.normalized()}
        , m_bias{Vector<Scalar, 3>::Zero()}
        , m_k_p{config.k_p}
        , m_k_i{config.k_i}
        , m_dt{config.dt}
        , m_gyro_buf{Vector<Scalar, 3>::Zero()}
        , m_health{cf_health::ok}
    {
        update_state_cache();
    }

    // Classifies without touching a member: the carried attitude and bias, then the sensors, then the step.
    ctrlpp::expected<void, cf_update_error> check_step(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt, const Vector<Scalar, 3>* mag = nullptr) const
    {
        if(!m_q.coeffs().allFinite() || !m_bias.allFinite())
            return ctrlpp::unexpected(cf_update_error::non_finite_state);
        if(!gyro.allFinite() || !accel.allFinite() || (mag != nullptr && !mag->allFinite()))
            return ctrlpp::unexpected(cf_update_error::non_finite_measurement);
        if(!std::isfinite(dt))
            return ctrlpp::unexpected(cf_update_error::non_finite_timestep);
        return {};
    }

    cf_update_error latch_health(cf_update_error fault)
    {
        if(fault == cf_update_error::non_finite_state)
            m_health = cf_health::non_finite_estimate;
        return fault;
    }

    // @cite mahony2008 -- Mahony, Hamel & Pflimlin, "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008, Eq. 48
    void integrate_gyro(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& e, Scalar dt)
    {
        m_bias -= m_k_i * e * dt;
        Vector<Scalar, 3> omega_c = gyro - m_bias + m_k_p * e;
        Vector<Scalar, 3> phi = dt * omega_c;
        m_q = (m_q * so3::exp(phi)).normalized();
        update_state_cache();
    }

    void update_state_cache() { m_state_cache << m_q.w(), m_q.x(), m_q.y(), m_q.z(), m_bias(0), m_bias(1), m_bias(2); }
};

static_assert(ObserverPolicy<complementary_filter<double>>);

}

#endif
