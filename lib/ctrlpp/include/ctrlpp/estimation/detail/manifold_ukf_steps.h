#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_MANIFOLD_UKF_STEPS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_MANIFOLD_UKF_STEPS_H

#include "ctrlpp/estimation/detail/geodesic_mean.h"

#include "ctrlpp/estimation/manifold_ukf_types.h"

#include "ctrlpp/types.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Cholesky>
#include <Eigen/Geometry>

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

// @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013
template <typename Scalar, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy>
class manifold_ukf_steps
{
public:
    static constexpr int ny = static_cast<int>(NY);
    static constexpr std::size_t num_sigma = Strategy::num_points;

    using state_vector_t = Vector<Scalar, 4>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, 3, 3>;
    using meas_cov_t = Matrix<Scalar, NY, NY>;

protected:
    Scalar m_tol;
    meas_cov_t m_R;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    Strategy m_strategy;
    Dynamics m_dynamics;
    Measurement m_measurement;
    std::int32_t m_max_iter;
    Eigen::Quaternion<Scalar> m_q;
    state_vector_t m_state_cache;
    output_vector_t m_innovation;
    manifold_ukf_health m_health;

    manifold_ukf_steps(Dynamics&& dynamics, Measurement&& measurement, manifold_ukf_config<Scalar, NY>&& config, Strategy&& strategy)
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

    template <typename Sigma>
    cov_matrix_t reduction_base(const Sigma& sigma) const
    {
        if(sigma.spd_repaired)
            return (sigma.factor * sigma.factor.transpose()).eval();
        return m_P;
    }

    void update_state_cache() { m_state_cache = so3::to_vec(m_q); }

    // @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013, Sec. 3.2.3
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

    // @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013, Sec. 3.2.3, Eq. 18
    detail::geodesic_mean_result<Scalar> compute_manifold_predicted_mean(const std::array<Eigen::Quaternion<Scalar>, num_sigma>& q_prop, const std::array<Scalar, num_sigma>& Wm) const
    {
        return detail::geodesic_mean_impl(q_prop, Wm, m_max_iter, m_tol);
    }

    // @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013, Sec. 3.2.3, Eq. 19
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

    // @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013, Sec. 3.2.4
    output_vector_t compute_predicted_measurement(const std::array<Scalar, num_sigma>& Wm, const std::array<output_vector_t, num_sigma>& z_sigma) const
    {
        output_vector_t z_pred = output_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_pred += Wm[i] * z_sigma[i];
        return z_pred;
    }

    // @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013, Sec. 3.2.4, Eqs. 36-37
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

}

#endif
