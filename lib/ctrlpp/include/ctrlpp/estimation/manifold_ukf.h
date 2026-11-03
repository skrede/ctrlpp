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

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/sigma_points/so3_sigma_points.h"

#include "ctrlpp/estimation/detail/geodesic_mean.h"
#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <Eigen/Cholesky>
#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy = so3_merwe_sigma_points<Scalar>>
    requires manifold_ukf_dynamics_model<Dynamics, Scalar> && manifold_ukf_measurement_model<Measurement, Scalar, NY> && manifold_sigma_point_strategy<Strategy, Scalar>
class manifold_ukf
{
    static_assert(NY > 0, "Output dimension NY must be positive");

    static constexpr int ny = static_cast<int>(NY);
    static constexpr std::size_t num_sigma = Strategy::num_points;

    struct validated_tag
    {
    };

public:
    using observer_tag = struct manifold_ukf_tag;
    using state_vector_t = Vector<Scalar, 4>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, 3, 3>;
    using meas_cov_t = Matrix<Scalar, NY, NY>;

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
        auto sigma = m_strategy.generate(m_q, m_P);
        note_health(sigma);
        auto q_prop = propagate_manifold_sigma_points(sigma.points, omega);
        auto mean_result = compute_manifold_predicted_mean(q_prop, sigma.Wm);
        auto q_new = mean_result.mean;
        m_P = compute_manifold_predicted_covariance(q_prop, q_new, sigma.Wc);
        m_q = q_new;
        update_state_cache();
        note_prediction(mean_result.converged);
    }

    // A refused step leaves every member bitwise unchanged, so the caller may retry with the next sample.
    ctrlpp::expected<void, manifold_ukf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = m_q;
        auto const previous_P = m_P;
        auto const previous_innovation = m_innovation;
        auto const previous_state = m_state_cache;
        auto sigma = m_strategy.generate(m_q, m_P);
        auto z_sigma = compute_measurement_sigma_points(sigma.points);
        auto z_pred = compute_predicted_measurement(sigma.Wm, z_sigma);
        auto [S, Pxz] = compute_innovation_and_cross_covariance(sigma, z_sigma, z_pred);
        auto K = compute_kalman_gain(Pxz, S);

        m_innovation = (z - z_pred).eval();
        apply_manifold_correction(K, S, reduction_base(sigma));
        update_state_cache();
        if(!m_q.coeffs().allFinite() || !m_P.allFinite() || !m_innovation.allFinite() || !m_state_cache.allFinite())
            return roll_back(previous_q, previous_P, previous_innovation, previous_state);
        note_health(sigma);
        note_definiteness();
        return {};
    }

    const state_vector_t& state() const { return m_state_cache; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

    const Eigen::Quaternion<Scalar>& attitude() const { return m_q; }

    manifold_ukf_health health() const { return m_health; }

private:
    Scalar m_tol;
    meas_cov_t m_R;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    Strategy m_strategy;
    Dynamics m_dynamics;
    Measurement m_measurement;
    std::size_t m_max_iter;
    Eigen::Quaternion<Scalar> m_q;
    state_vector_t m_state_cache;
    output_vector_t m_innovation;
    manifold_ukf_health m_health;

    manifold_ukf(validated_tag, Dynamics dynamics, Measurement measurement, manifold_ukf_config<Scalar, NY> config, Strategy strategy)
        : m_tol{config.geodesic_mean_tol}
        , m_R{std::move(config.R)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_strategy{std::move(strategy)}
        , m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_max_iter{config.geodesic_mean_max_iter}
        , m_q{config.q0.normalized()}
        , m_state_cache{so3::to_vec(m_q)}
        , m_innovation{output_vector_t::Zero()}
        , m_health{manifold_ukf_health::ok}
    {
    }

    ctrlpp::expected<void, manifold_ukf_update_error> roll_back(
            const Eigen::Quaternion<Scalar>& q, const cov_matrix_t& P, const output_vector_t& innovation, const state_vector_t& state)
    {
        m_q = q;
        m_P = P;
        m_innovation = innovation;
        m_state_cache = state;
        return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_result);
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, manifold_ukf_update_error> check_step(const output_vector_t& z) const
    {
        if(!m_q.coeffs().allFinite())
            return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_state);
        if(!m_P.allFinite())
            return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(manifold_ukf_update_error::non_finite_measurement);
        return {};
    }

    manifold_ukf_update_error latch_health(manifold_ukf_update_error fault)
    {
        if(fault != manifold_ukf_update_error::non_finite_measurement)
            m_health = manifold_ukf_health::non_finite_estimate;
        return fault;
    }

    template <typename Sigma>
    void note_health(const Sigma& sigma)
    {
        if(sigma.spd_repaired && m_health == manifold_ukf_health::ok)
            m_health = manifold_ukf_health::covariance_repaired;
    }

    void note_prediction(bool mean_converged)
    {
        if(!m_q.coeffs().allFinite() || !m_P.allFinite())
            m_health = manifold_ukf_health::non_finite_estimate;
        else if(!mean_converged && m_health < manifold_ukf_health::mean_not_converged)
            m_health = manifold_ukf_health::mean_not_converged;
        else
            note_definiteness();
    }

    // Judged by the unpivoted Cholesky factorization the next sigma-point generation attempts first on this
    // same matrix, so the status never disagrees with the repair that generation makes.
    void note_definiteness()
    {
        if(m_health == manifold_ukf_health::ok && Eigen::LLT<cov_matrix_t>{m_P}.info() != Eigen::Success)
            m_health = manifold_ukf_health::covariance_repaired;
    }

    template <typename Sigma>
    cov_matrix_t reduction_base(const Sigma& sigma) const
    {
        if(sigma.spd_repaired)
            return (sigma.factor * sigma.factor.transpose()).eval();
        return m_P;
    }

    void update_state_cache() { m_state_cache = so3::to_vec(m_q); }

    // @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013, Sec. 3.2
    std::array<Eigen::Quaternion<Scalar>, num_sigma> propagate_manifold_sigma_points(const std::array<Eigen::Quaternion<Scalar>, num_sigma>& points, const input_vector_t& omega) const
    {
        std::array<Eigen::Quaternion<Scalar>, num_sigma> q_prop;
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            q_prop[i] = m_dynamics(points[i], omega);
            q_prop[i].normalize();
        }
        return q_prop;
    }

    // @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013, Alg. 1
    detail::geodesic_mean_result<Scalar> compute_manifold_predicted_mean(const std::array<Eigen::Quaternion<Scalar>, num_sigma>& q_prop, const std::array<Scalar, num_sigma>& Wm) const
    {
        return detail::geodesic_mean_impl(q_prop, Wm, m_max_iter, m_tol);
    }

    // @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013, Eq. 10
    cov_matrix_t compute_manifold_predicted_covariance(
            const std::array<Eigen::Quaternion<Scalar>, num_sigma>& q_prop, const Eigen::Quaternion<Scalar>& q_mean, const std::array<Scalar, num_sigma>& Wc) const
    {
        cov_matrix_t P_pred = cov_matrix_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            Eigen::Quaternion<Scalar> qi = q_prop[i];
            if(qi.dot(q_mean) < Scalar{0})
            {
                qi.w() = -qi.w();
                qi.x() = -qi.x();
                qi.y() = -qi.y();
                qi.z() = -qi.z();
            }
            Vector<Scalar, 3> e = so3::log(q_mean.conjugate() * qi);
            P_pred += Wc[i] * e * e.transpose();
        }
        P_pred += m_Q;
        return detail::symmetrize(P_pred);
    }

    std::array<output_vector_t, num_sigma> compute_measurement_sigma_points(const std::array<Eigen::Quaternion<Scalar>, num_sigma>& points) const
    {
        std::array<output_vector_t, num_sigma> z_sigma;
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_sigma[i] = m_measurement(points[i]);
        return z_sigma;
    }

    // @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013
    output_vector_t compute_predicted_measurement(const std::array<Scalar, num_sigma>& Wm, const std::array<output_vector_t, num_sigma>& z_sigma) const
    {
        output_vector_t z_pred = output_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_pred += Wm[i] * z_sigma[i];
        return z_pred;
    }

    // @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013, Eq. 12-13
    std::pair<meas_cov_t, Eigen::Matrix<Scalar, 3, ny>> compute_innovation_and_cross_covariance(
            const auto& sigma, const std::array<output_vector_t, num_sigma>& z_sigma, const output_vector_t& z_pred) const
    {
        meas_cov_t S = meas_cov_t::Zero();
        Eigen::Matrix<Scalar, 3, ny> Pxz = Eigen::Matrix<Scalar, 3, ny>::Zero();

        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto dz = (z_sigma[i] - z_pred).eval();
            S += sigma.Wc[i] * dz * dz.transpose();

            Eigen::Quaternion<Scalar> qi = sigma.points[i];
            if(qi.dot(m_q) < Scalar{0})
            {
                qi.w() = -qi.w();
                qi.x() = -qi.x();
                qi.y() = -qi.y();
                qi.z() = -qi.z();
            }
            Vector<Scalar, 3> dx = so3::log(m_q.conjugate() * qi);
            Pxz += sigma.Wc[i] * dx * dz.transpose();
        }
        S += m_R;

        return {S, Pxz};
    }

    // @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
    Eigen::Matrix<Scalar, 3, ny> compute_kalman_gain(const Eigen::Matrix<Scalar, 3, ny>& Pxz, const meas_cov_t& S) const
    {
        if constexpr(ny == 1)
            return detail::one_output_gain_transpose(S(0, 0), Pxz.transpose()).transpose().eval();
        else
        {
            Eigen::Matrix<Scalar, ny, 3> KT = S.transpose().colPivHouseholderQr().solve(Pxz.transpose());
            return KT.transpose().eval();
        }
    }

    // The reduction starts from the tangent covariance the sigma set represents, the repaired one after a
    // repair, since K and S describe that one.
    // @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
    void apply_manifold_correction(const Eigen::Matrix<Scalar, 3, ny>& K, const meas_cov_t& S, const cov_matrix_t& base)
    {
        Vector<Scalar, 3> delta_phi = K * m_innovation;
        m_q = (m_q * so3::exp(delta_phi)).normalized();

        // The reset Jacobian G = I - 0.5 skew(delta_phi), as in the MEKF reset, carries the covariance into the
        // corrected frame; without it the covariance stays anchored to the pre-correction frame.
        cov_matrix_t G = cov_matrix_t::Identity() - Scalar{0.5} * so3::skew(delta_phi);
        m_P = detail::symmetrize((G * (base - K * S * K.transpose()) * G.transpose()).eval());
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
