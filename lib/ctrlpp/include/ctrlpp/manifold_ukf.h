#ifndef HPP_GUARD_CTRLPP_MANIFOLD_UKF_H
#define HPP_GUARD_CTRLPP_MANIFOLD_UKF_H

#include "ctrlpp/so3.h"
#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"
#include "ctrlpp/sigma_points/so3_sigma_points.h"

#include <Eigen/Geometry>

#include <array>
#include <cstddef>
#include <utility>

namespace ctrlpp {

template <typename M, typename Scalar, std::size_t NY>
concept manifold_ukf_measurement_model = requires(const M &m, const Eigen::Quaternion<Scalar> &q)
{
    { m(q) } -> std::convertible_to<Vector<Scalar, NY>>;
};

template <typename D, typename Scalar>
concept manifold_ukf_dynamics_model = requires(const D &d, const Eigen::Quaternion<Scalar> &q, const Vector<Scalar, 3> &u)
{
    { d(q, u) } -> std::convertible_to<Eigen::Quaternion<Scalar>>;
};

template <typename Scalar, std::size_t NY>
struct manifold_ukf_config
{
    Matrix<Scalar, 3, 3> Q{Matrix<Scalar, 3, 3>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
    Matrix<Scalar, 3, 3> P0{Matrix<Scalar, 3, 3>::Identity()};
    Scalar dt{Scalar{0.01}};
    std::size_t geodesic_mean_max_iter{30};
    Scalar geodesic_mean_tol{Scalar{1e-9}};
};

namespace detail {

template <typename Scalar, std::size_t NP>
Eigen::Quaternion<Scalar> geodesic_mean_impl(const std::array<Eigen::Quaternion<Scalar>, NP> &qs, const std::array<Scalar, NP> &Wm, std::size_t max_iter, Scalar tol)
{
    auto q_mean = qs[0];
    for(std::size_t iter = 0; iter < max_iter; ++iter)
    {
        Vector<Scalar, 3> eps = Vector<Scalar, 3>::Zero();
        for(std::size_t i = 0; i < NP; ++i)
        {
            Eigen::Quaternion<Scalar> qi = qs[i];
            if(qi.dot(q_mean) < Scalar{0})
            {
                qi.w() = -qi.w();
                qi.x() = -qi.x();
                qi.y() = -qi.y();
                qi.z() = -qi.z();
            }
            eps += Wm[i] * so3::log(q_mean.conjugate() * qi);
        }
        if(eps.norm() < tol)
            break;
        q_mean = (q_mean * so3::exp(eps)).normalized();
    }
    return q_mean;
}

}

template <typename Scalar, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy = so3_merwe_sigma_points<Scalar>>
    requires manifold_ukf_dynamics_model<Dynamics, Scalar> && manifold_ukf_measurement_model<Measurement, Scalar, NY> && manifold_sigma_point_strategy<Strategy, Scalar>
class manifold_ukf
{
    static constexpr int ny = static_cast<int>(NY);
    static constexpr std::size_t num_sigma = Strategy::num_points;

public:
    using observer_tag = struct manifold_ukf_tag;
    using state_vector_t = Vector<Scalar, 4>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, 3, 3>;
    using meas_cov_t = Matrix<Scalar, NY, NY>;

    manifold_ukf(Dynamics dynamics, Measurement measurement, manifold_ukf_config<Scalar, NY> config, Strategy strategy = Strategy{})
        : m_tol{config.geodesic_mean_tol}
        , m_R{std::move(config.R)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_strategy{std::move(strategy)}
        , m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_max_iter{config.geodesic_mean_max_iter}
        , m_state_cache{so3::to_vec(m_q)}
        , m_innovation{output_vector_t::Zero()}
        , m_q{config.q0.normalized()}
    {
    }

    void predict(const input_vector_t &omega)
    {
        auto sigma = m_strategy.generate(m_q, m_P);

        std::array<Eigen::Quaternion<Scalar>, num_sigma> q_prop;
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            q_prop[i] = dynamics_(sigma.points[i], omega);
            q_prop[i].normalize();
        }

        auto q_new = detail::geodesic_mean_impl(q_prop, sigma.Wm, m_max_iter, m_tol);

        cov_matrix_t P_pred = cov_matrix_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            Eigen::Quaternion<Scalar> qi = q_prop[i];
            if(qi.dot(q_new) < Scalar{0})
            {
                qi.w() = -qi.w();
                qi.x() = -qi.x();
                qi.y() = -qi.y();
                qi.z() = -qi.z();
            }
            Vector<Scalar, 3> e = so3::log(q_new.conjugate() * qi);
            P_pred += sigma.Wc[i] * e * e.transpose();
        }
        P_pred += m_Q;

        m_q = q_new;
        m_P = (Scalar{0.5} * (P_pred + P_pred.transpose())).eval();
        update_state_cache();
    }

    void update(const output_vector_t &z)
    {
        auto sigma = m_strategy.generate(m_q, m_P);

        std::array<output_vector_t, num_sigma> z_sigma;
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_sigma[i] = measurement_(sigma.points[i]);

        output_vector_t z_pred = output_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_pred += sigma.Wm[i] * z_sigma[i];

        meas_cov_t S = meas_cov_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto dz = (z_sigma[i] - z_pred).eval();
            S += sigma.Wc[i] * dz * dz.transpose();
        }
        S += m_R;

        Eigen::Matrix<Scalar, 3, ny> Pxz = Eigen::Matrix<Scalar, 3, ny>::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            Eigen::Quaternion<Scalar> qi = sigma.points[i];
            if(qi.dot(m_q) < Scalar{0})
            {
                qi.w() = -qi.w();
                qi.x() = -qi.x();
                qi.y() = -qi.y();
                qi.z() = -qi.z();
            }
            Vector<Scalar, 3> dx = so3::log(m_q.conjugate() * qi);
            auto dz = (z_sigma[i] - z_pred).eval();
            Pxz += sigma.Wc[i] * dx * dz.transpose();
        }

        Eigen::Matrix<Scalar, ny, 3> KT = S.transpose().colPivHouseholderQr().solve(Pxz.transpose());
        Eigen::Matrix<Scalar, 3, ny> K = KT.transpose().eval();

        m_innovation = (z - z_pred).eval();

        Vector<Scalar, 3> delta_phi = K * m_innovation;
        m_q = (m_q * so3::exp(delta_phi)).normalized();

        cov_matrix_t P_new = (m_P - K * S * K.transpose()).eval();
        m_P = (Scalar{0.5} * (P_new + P_new.transpose())).eval();
        update_state_cache();
    }

    const state_vector_t &state() const
    {
        return m_state_cache;
    }

    const cov_matrix_t &covariance() const
    {
        return m_P;
    }

    const output_vector_t &innovation() const
    {
        return m_innovation;
    }

    const Eigen::Quaternion<Scalar> &attitude() const
    {
        return m_q;
    }

private:
    Scalar m_tol;
    meas_cov_t m_R;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    Strategy m_strategy;
    Dynamics m_dynamics;
    Measurement m_measurement;
    std::size_t m_max_iter;
    state_vector_t m_state_cache;
    output_vector_t m_innovation;
    Eigen::Quaternion<Scalar> m_q;

    void update_state_cache()
    {
        m_state_cache = so3::to_vec(m_q);
    }

};

template <typename Dynamics, typename Measurement, typename Scalar, std::size_t NY>
manifold_ukf(Dynamics, Measurement, manifold_ukf_config<Scalar, NY>) -> manifold_ukf<Scalar, NY, Dynamics, Measurement, so3_merwe_sigma_points<Scalar>>;

namespace detail {

struct mukf_sa_dynamics
{
    auto operator()(const Eigen::Quaternion<double> &q, const Vector<double, 3> &) const -> Eigen::Quaternion<double> { return q; }
};

struct mukf_sa_measurement
{
    auto operator()(const Eigen::Quaternion<double> &) const -> Vector<double, 3> { return Vector<double, 3>::Zero(); }
};

}

static_assert(ObserverPolicy<manifold_ukf<double, 3, detail::mukf_sa_dynamics, detail::mukf_sa_measurement>>);
static_assert(CovarianceObserver<manifold_ukf<double, 3, detail::mukf_sa_dynamics, detail::mukf_sa_measurement>>);

}

#endif
