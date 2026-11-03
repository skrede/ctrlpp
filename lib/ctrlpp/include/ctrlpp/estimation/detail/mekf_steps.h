#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_MEKF_STEPS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_MEKF_STEPS_H

#include "ctrlpp/estimation/mekf_types.h"

#include "ctrlpp/types.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/detail/covariance_ops.h"
#include "ctrlpp/detail/numerical_mekf_diff.h"

#include <Eigen/Geometry>

#include <cstddef>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NB, std::size_t NY, typename Measurement>
class mekf_steps
{
public:
    static constexpr std::size_t NE = 3 + NB;
    static constexpr int ne = static_cast<int>(NE);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr int nb = static_cast<int>(NB);

    using state_vector_t = Vector<Scalar, 4 + NB>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NE, NE>;
    using meas_cov_t = Matrix<Scalar, NY, NY>;

protected:
    Measurement m_measurement;
    Eigen::Quaternion<Scalar> m_q;
    Vector<Scalar, NB> m_b;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    meas_cov_t m_R;
    Scalar m_eps;
    Scalar m_dt;
    output_vector_t m_innovation;
    state_vector_t m_state_cache;
    mekf_health m_health;

    mekf_steps(Measurement&& measurement, mekf_config<Scalar, NB, NY>&& config)
        : m_measurement{std::move(measurement)}
        , m_q{config.q0.normalized()}
        , m_b{std::move(config.b0)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_eps{config.numerical_eps}
        , m_dt{config.dt}
        , m_innovation{output_vector_t::Zero()}
        , m_state_cache{}
        , m_health{mekf_health::ok}
    {
        update_state_cache();
    }

    // @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
    // @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
    void predict_impl(const input_vector_t& omega, Scalar dt)
    {
        Vector<Scalar, 3> omega_corr = omega - m_b.template head<3>();
        Vector<Scalar, 3> omega_dt = (omega_corr * dt).eval();
        m_q = (m_q * so3::exp(omega_dt)).normalized();

        Eigen::Matrix<Scalar, 3, 3> C = so3::exp(omega_dt).toRotationMatrix();

        // The error is right (body-frame) multiplicative, q = q_nominal * exp(delta), so it propagates
        // through exp(-omega_dt) * exp(delta) * exp(omega_dt) = exp(C^T delta): the attitude block is C^T.
        cov_matrix_t F = cov_matrix_t::Identity();
        F.template block<3, 3>(0, 0) = C.transpose();
        F.template block<3, nb>(0, 3) = -Eigen::Matrix<Scalar, 3, nb>::Identity() * dt;

        m_P = detail::symmetrize((F * m_P * F.transpose() + m_Q).eval());
        update_state_cache();
    }

    // @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003, Eq. 34
    Matrix<Scalar, NY, NE> compute_measurement_jacobian() const
    {
        if constexpr(differentiable_mekf_measurement<Measurement, Scalar, NB, NY>)
            return m_measurement.jacobian(m_q, m_b);
        else
            return detail::numerical_mekf_jacobian<Scalar, NB, NY>(m_measurement, m_q, m_b, m_eps);
    }

    meas_cov_t compute_innovation_covariance(const Matrix<Scalar, NY, NE>& H) const { return (H * m_P * H.transpose() + m_R).eval(); }

    // @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
    Eigen::Matrix<Scalar, ne, ny> compute_kalman_gain(const Matrix<Scalar, NY, NE>& H, const meas_cov_t& S) const
    {
        if constexpr(ny == 1)
            return detail::one_output_gain_transpose(S(0, 0), H * m_P).transpose().eval();
        else
        {
            Eigen::Matrix<Scalar, ny, ne> KT = S.transpose().colPivHouseholderQr().solve(H * m_P);
            return KT.transpose().eval();
        }
    }

    // @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003, Eq. 46
    Vector<Scalar, NE> apply_multiplicative_correction(const Eigen::Matrix<Scalar, ne, ny>& K)
    {
        Vector<Scalar, NE> delta_xi = K * m_innovation;
        Vector<Scalar, 3> delta_att = delta_xi.template head<3>().eval();
        m_q = (m_q * so3::exp(delta_att)).normalized();
        m_b += delta_xi.template tail<NB>();
        return delta_xi;
    }

    // Joseph form, then the reset through the frame-change Jacobian G that the multiplicative correction
    // makes mandatory.
    // @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003, Eq. 68
    void update_covariance(const Eigen::Matrix<Scalar, ne, ny>& K, const Matrix<Scalar, NY, NE>& H, const Vector<Scalar, NE>& delta_xi)
    {
        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        m_P = (IKH * m_P * IKH.transpose() + K * m_R * K.transpose()).eval();

        Vector<Scalar, 3> delta_att = delta_xi.template head<3>();
        cov_matrix_t G = cov_matrix_t::Identity();
        G.template block<3, 3>(0, 0) -= Scalar{0.5} * so3::skew(delta_att);
        m_P = detail::symmetrize((G * m_P * G.transpose()).eval());
    }

    void update_state_cache()
    {
        m_state_cache(0) = m_q.w();
        m_state_cache(1) = m_q.x();
        m_state_cache(2) = m_q.y();
        m_state_cache(3) = m_q.z();
        for(std::size_t i = 0; i < NB; ++i)
        {
            m_state_cache(static_cast<Eigen::Index>(4 + i)) = m_b(static_cast<Eigen::Index>(i));
        }
    }
};

}

#endif
