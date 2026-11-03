#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_EKF_STEPS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_EKF_STEPS_H

#include "ctrlpp/estimation/ekf_types.h"

#include "ctrlpp/types.h"

#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"

#include "ctrlpp/detail/covariance_ops.h"
#include "ctrlpp/detail/numerical_diff.h"

#include <cstddef>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement>
class ekf_steps
{
public:
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);

    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;

protected:
    Scalar m_eps;
    Scalar m_nis;
    Dynamics m_dynamics;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    state_vector_t m_x;
    meas_cov_matrix_t m_R;
    Measurement m_measurement;
    output_vector_t m_innovation;
    ekf_health m_health;

    ekf_steps(Dynamics&& dynamics, Measurement&& measurement, ekf_config<Scalar, NX, NU, NY>&& config)
        : m_eps{config.numerical_eps}
        , m_nis{0}
        , m_dynamics{std::move(dynamics)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_x{std::move(config.x0)}
        , m_R{std::move(config.R)}
        , m_measurement{std::move(measurement)}
        , m_innovation{output_vector_t::Zero()}
        , m_health{ekf_health::ok}
    {
    }

    void propagate_state(const state_vector_t& x_prev, const input_vector_t& u) { m_x = m_dynamics(x_prev, u); }

    // @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.1
    void propagate_covariance(const state_vector_t& x_prev, const input_vector_t& u)
    {
        Matrix<Scalar, NX, NX> F;
        if constexpr(differentiable_dynamics<Dynamics, Scalar, NX, NU>)
            F = m_dynamics.jacobian_x(x_prev, u);
        else
            F = detail::numerical_jacobian_x<Scalar, NX, NU>(m_dynamics, x_prev, u, m_eps);

        m_P = detail::symmetrize((F * m_P * F.transpose() + m_Q).eval());
    }

    Matrix<Scalar, NY, NX> compute_measurement_jacobian() const
    {
        if constexpr(differentiable_measurement<Measurement, Scalar, NX, NY>)
            return m_measurement.jacobian(m_x);
        else
            return detail::numerical_jacobian_h<Scalar, NX, NY>(m_measurement, m_x, m_eps);
    }

    // @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.2
    meas_cov_matrix_t compute_innovation_covariance(const Matrix<Scalar, NY, NX>& H) const { return (H * m_P * H.transpose() + m_R).eval(); }

    // @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.2
    Eigen::Matrix<Scalar, nx, ny> compute_kalman_gain(const Matrix<Scalar, NY, NX>& H, const meas_cov_matrix_t& S) const
    {
        if constexpr(ny == 1)
            return detail::one_output_gain_transpose(S(0, 0), H * m_P).transpose().eval();
        else
        {
            Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(H * m_P);
            return KT_solved.transpose().eval();
        }
    }

    void apply_state_correction(const Eigen::Matrix<Scalar, nx, ny>& K, const output_vector_t& /*z*/, const output_vector_t& /*z_pred*/) { m_x = (m_x + K * m_innovation).eval(); }

    // Joseph form, symmetrized.
    // @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.3
    void update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K, const Matrix<Scalar, NY, NX>& H)
    {
        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        m_P = detail::symmetrize((IKH * m_P * IKH.transpose() + K * m_R * K.transpose()).eval());
    }
};

}

#endif
