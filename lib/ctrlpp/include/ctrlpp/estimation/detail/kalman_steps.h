#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_KALMAN_STEPS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_KALMAN_STEPS_H

#include "ctrlpp/estimation/kalman_types.h"

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <cstddef>
#include <utility>

namespace ctrlpp::detail
{

// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class kalman_steps
{
public:
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

    using state_vector_t = Eigen::Matrix<Scalar, nx, 1>;
    using input_vector_t = Eigen::Matrix<Scalar, nu, 1>;
    using output_vector_t = Eigen::Matrix<Scalar, ny, 1>;
    using cov_matrix_t = Eigen::Matrix<Scalar, nx, nx>;
    using meas_cov_matrix_t = Eigen::Matrix<Scalar, ny, ny>;
    using system_t = discrete_state_space<Scalar, NX, NU, NY>;

protected:
    Scalar m_nis_value;
    system_t m_sys;
    cov_matrix_t m_Q;
    cov_matrix_t m_P;
    cov_matrix_t m_P_post_prev;
    state_vector_t m_x;
    meas_cov_matrix_t m_R;
    output_vector_t m_innovation;
    kalman_health m_health;

    kalman_steps(system_t&& sys, kalman_config<Scalar, NX, NU, NY>&& config)
        : m_nis_value{0}
        , m_sys{std::move(sys)}
        , m_Q{std::move(config.Q)}
        , m_P{std::move(config.P0)}
        , m_P_post_prev{m_P}
        , m_x{std::move(config.x0)}
        , m_R{std::move(config.R)}
        , m_innovation{output_vector_t::Zero()}
        , m_health{kalman_health::ok}
    {
    }

    void propagate_state(const input_vector_t& u) { m_x = (m_sys.A * m_x + m_sys.B * u).eval(); }

    void propagate_covariance() { m_P = (m_sys.A * m_P * m_sys.A.transpose() + m_Q).eval(); }

    void compute_innovation(const output_vector_t& z) { m_innovation = (z - m_sys.C * m_x).eval(); }

    meas_cov_matrix_t compute_innovation_covariance() const { return (m_sys.C * m_P * m_sys.C.transpose() + m_R).eval(); }

    Eigen::Matrix<Scalar, nx, ny> compute_kalman_gain(const meas_cov_matrix_t& S) const
    {
        Eigen::Matrix<Scalar, ny, nx> CP = m_sys.C * m_P;
        if constexpr(ny == 1)
            return detail::one_output_gain_transpose(S(0, 0), CP).transpose().eval();
        else
        {
            Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(CP);
            return KT_solved.transpose().eval();
        }
    }

    void apply_state_correction(const Eigen::Matrix<Scalar, nx, ny>& K) { m_x = (m_x + K * m_innovation).eval(); }

    // Joseph form, symmetrized.
    void update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K)
    {
        cov_matrix_t IKC = cov_matrix_t::Identity() - K * m_sys.C;
        m_P = detail::symmetrize((IKC * m_P * IKC.transpose() + K * m_R * K.transpose()).eval());
    }

    void compute_nis(const meas_cov_matrix_t& S)
    {
        output_vector_t Sinv_z = S.colPivHouseholderQr().solve(m_innovation).eval();
        m_nis_value = (m_innovation.transpose() * Sinv_z)(0, 0);
    }
};

}

#endif
