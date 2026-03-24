#ifndef HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H
#define HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H

/// @brief Linear discrete-time Kalman filter with Joseph-form covariance update.

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/detail/covariance_ops.h"

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

    /// @brief Predict state and covariance one step forward.
    void predict(const input_vector_t& u)
    {
        m_P_post_prev = m_P;
        propagate_state(u);
        propagate_covariance();
    }

    /// @brief Update state and covariance with measurement.
    void update(const output_vector_t& z)
    {
        compute_innovation(z);
        auto S = compute_innovation_covariance();
        auto K = compute_kalman_gain(S);

        apply_state_correction(K);
        update_covariance(K);
        compute_nees(S);
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
    /// @brief Propagate state: x = A*x + B*u.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    void propagate_state(const input_vector_t& u)
    {
        m_x = (m_sys.A * m_x + m_sys.B * u).eval();
    }

    /// @brief Propagate covariance: P = A*P*A^T + Q.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    void propagate_covariance()
    {
        m_P = (m_sys.A * m_P * m_sys.A.transpose() + m_Q).eval();
    }

    /// @brief Compute innovation: y = z - C*x.
    void compute_innovation(const output_vector_t& z)
    {
        m_innovation = (z - m_sys.C * m_x).eval();
    }

    /// @brief Compute innovation covariance: S = C*P*C^T + R.
    [[nodiscard]] auto compute_innovation_covariance() const -> meas_cov_matrix_t
    {
        return (m_sys.C * m_P * m_sys.C.transpose() + m_R).eval();
    }

    /// @brief Compute Kalman gain: K = P*C^T*S^{-1} via column-pivoting QR solve.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    [[nodiscard]] auto compute_kalman_gain(const meas_cov_matrix_t& S) const -> Eigen::Matrix<Scalar, nx, ny>
    {
        Eigen::Matrix<Scalar, ny, nx> CP = m_sys.C * m_P;
        Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(CP);
        return KT_solved.transpose().eval();
    }

    /// @brief Apply state correction: x += K * innovation.
    void apply_state_correction(const Eigen::Matrix<Scalar, nx, ny>& K)
    {
        m_x = (m_x + K * m_innovation).eval();
    }

    /// @brief Joseph-form covariance update: P = (I-KC)*P*(I-KC)^T + K*R*K^T, symmetrized.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    void update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K)
    {
        cov_matrix_t IKC = cov_matrix_t::Identity() - K * m_sys.C;
        m_P = detail::symmetrize((IKC * m_P * IKC.transpose() + K * m_R * K.transpose()).eval());
    }

    /// @brief Compute Normalized Estimation Error Squared: y^T * S^{-1} * y.
    void compute_nees(const meas_cov_matrix_t& S)
    {
        output_vector_t Sinv_z = S.colPivHouseholderQr().solve(m_innovation).eval();
        m_nees_value = (m_innovation.transpose() * Sinv_z)(0, 0);
    }

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

}

#endif
