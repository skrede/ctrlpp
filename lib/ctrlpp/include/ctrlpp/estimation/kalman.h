#ifndef HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H
#define HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H

/// @brief Linear discrete-time Kalman filter with Joseph-form covariance update.
///
/// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/estimation/observer_policy.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct kalman_config
{
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class kalman_filter
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct kalman_tag;
    using state_vector_t = Eigen::Matrix<Scalar, nx, 1>;
    using input_vector_t = Eigen::Matrix<Scalar, nu, 1>;
    using output_vector_t = Eigen::Matrix<Scalar, ny, 1>;
    using cov_matrix_t = Eigen::Matrix<Scalar, nx, nx>;
    using meas_cov_matrix_t = Eigen::Matrix<Scalar, ny, ny>;
    using system_t = discrete_state_space<Scalar, NX, NU, NY>;

    kalman_filter(system_t sys, kalman_config<Scalar, NX, NU, NY> config)
        : m_sys{std::move(sys)}, m_Q{std::move(config.Q)}, m_P{std::move(config.P0)}, m_P_post_prev{m_P}, m_x{std::move(config.x0)}, m_R{std::move(config.R)}, m_innovation{output_vector_t::Zero()}
    {
    }

    void predict(const input_vector_t& u)
    {
        // Save current P (post-update from previous cycle) for steady-state comparison
        m_P_post_prev = m_P;
        m_x = (m_sys.A * m_x + m_sys.B * u).eval();
        m_P = (m_sys.A * m_P * m_sys.A.transpose() + m_Q).eval();
    }

    void update(const output_vector_t& z)
    {
        // Innovation
        m_innovation = (z - m_sys.C * m_x).eval();

        // Innovation covariance S = C P C^T + R
        meas_cov_matrix_t S = (m_sys.C * m_P * m_sys.C.transpose() + m_R).eval();

        // Kalman gain: K = P C^T S^{-1}
        // Solve via: S^T K^T = C P
        Eigen::Matrix<Scalar, ny, nx> CP = m_sys.C * m_P;
        Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(CP);
        Eigen::Matrix<Scalar, nx, ny> K = KT_solved.transpose().eval();

        // State update
        m_x = (m_x + K * m_innovation).eval();

        // Joseph form covariance update: P = (I - KC) P (I - KC)^T + K R K^T
        cov_matrix_t IKC = cov_matrix_t::Identity() - K * m_sys.C;
        cov_matrix_t P_new = (IKC * m_P * IKC.transpose() + K * m_R * K.transpose()).eval();

        // Symmetrize to combat floating point drift
        m_P = (Scalar{0.5} * (P_new + P_new.transpose())).eval();

        // NEES: z^T S^{-1} z
        output_vector_t Sinv_z = S.colPivHouseholderQr().solve(m_innovation).eval();
        m_nees_value = (m_innovation.transpose() * Sinv_z)(0, 0);
    }

    const state_vector_t& state() const { return m_x; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

    Scalar nees() const { return m_nees_value; }

    bool is_steady_state(Scalar tol = Scalar{1e-10}) const
    {
        Scalar p_norm = m_P.norm();
        if(p_norm < std::numeric_limits<Scalar>::epsilon())
            return true;
        return (m_P - m_P_post_prev).norm() / p_norm < tol;
    }

    void reset_covariance(const cov_matrix_t& P0)
    {
        m_P = P0;
        m_P_post_prev = P0;
    }

    void set_model(system_t sys) { m_sys = std::move(sys); }

    void set_noise(cov_matrix_t Q, meas_cov_matrix_t R)
    {
        m_Q = std::move(Q);
        m_R = std::move(R);
    }

private:
    Scalar m_nees_value{0};
    system_t m_sys;
    cov_matrix_t m_Q;
    cov_matrix_t m_P;
    cov_matrix_t m_P_post_prev;
    state_vector_t m_x;
    meas_cov_matrix_t m_R;
    output_vector_t m_innovation;
};

static_assert(ObserverPolicy<kalman_filter<double, 2, 1, 1>>);
static_assert(CovarianceObserver<kalman_filter<double, 2, 1, 1>>);

} // namespace ctrlpp

#endif
