#ifndef HPP_GUARD_CTRLPP_KALMAN_H
#define HPP_GUARD_CTRLPP_KALMAN_H

#include "ctrlpp/types.h"
#include "ctrlpp/state_space.h"
#include "ctrlpp/observer_policy.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class KalmanFilter {
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
    using system_t = DiscreteStateSpace<Scalar, NX, NU, NY>;

    KalmanFilter(system_t sys,
                 cov_matrix_t Q,
                 meas_cov_matrix_t R,
                 state_vector_t x0,
                 cov_matrix_t P0)
        : sys_{std::move(sys)}
        , x_{std::move(x0)}
        , P_{std::move(P0)}
        , P_post_prev_{P_}
        , Q_{std::move(Q)}
        , R_{std::move(R)}
        , innovation_{output_vector_t::Zero()}
    {
    }

    void predict(const input_vector_t& u)
    {
        // Save current P (post-update from previous cycle) for steady-state comparison
        P_post_prev_ = P_;
        x_ = (sys_.A * x_ + sys_.B * u).eval();
        P_ = (sys_.A * P_ * sys_.A.transpose() + Q_).eval();
    }

    void update(const output_vector_t& z)
    {
        // Innovation
        innovation_ = (z - sys_.C * x_).eval();

        // Innovation covariance S = C P C^T + R
        meas_cov_matrix_t S = (sys_.C * P_ * sys_.C.transpose() + R_).eval();

        // Kalman gain: K = P C^T S^{-1}
        // Solve via: S^T K^T = C P
        Eigen::Matrix<Scalar, ny, nx> CP = sys_.C * P_;
        Eigen::Matrix<Scalar, ny, nx> KT_solved =
            S.transpose().colPivHouseholderQr().solve(CP);
        Eigen::Matrix<Scalar, nx, ny> K = KT_solved.transpose().eval();

        // State update
        x_ = (x_ + K * innovation_).eval();

        // Joseph form covariance update: P = (I - KC) P (I - KC)^T + K R K^T
        cov_matrix_t IKC = cov_matrix_t::Identity() - K * sys_.C;
        cov_matrix_t P_new = (IKC * P_ * IKC.transpose() + K * R_ * K.transpose()).eval();

        // Symmetrize to combat floating point drift
        P_ = (Scalar{0.5} * (P_new + P_new.transpose())).eval();

        // NEES: z^T S^{-1} z
        output_vector_t Sinv_z = S.colPivHouseholderQr().solve(innovation_).eval();
        nees_value_ = (innovation_.transpose() * Sinv_z)(0, 0);
    }

    [[nodiscard]] auto state() const -> const state_vector_t& { return x_; }
    [[nodiscard]] auto covariance() const -> const cov_matrix_t& { return P_; }
    [[nodiscard]] auto innovation() const -> const output_vector_t& { return innovation_; }
    [[nodiscard]] auto nees() const -> Scalar { return nees_value_; }

    [[nodiscard]] auto is_steady_state(Scalar tol = Scalar{1e-10}) const -> bool
    {
        Scalar p_norm = P_.norm();
        if (p_norm < std::numeric_limits<Scalar>::epsilon())
            return true;
        return (P_ - P_post_prev_).norm() / p_norm < tol;
    }

    void reset_covariance(const cov_matrix_t& P0)
    {
        P_ = P0;
        P_post_prev_ = P0;
    }

    void set_model(system_t sys) { sys_ = std::move(sys); }

    void set_noise(cov_matrix_t Q, meas_cov_matrix_t R)
    {
        Q_ = std::move(Q);
        R_ = std::move(R);
    }

private:
    system_t sys_;
    state_vector_t x_;
    cov_matrix_t P_;
    cov_matrix_t P_post_prev_;
    cov_matrix_t Q_;
    meas_cov_matrix_t R_;
    output_vector_t innovation_;
    Scalar nees_value_{0};
};

static_assert(ObserverPolicy<KalmanFilter<double, 2, 1, 1>>);
static_assert(CovarianceObserver<KalmanFilter<double, 2, 1, 1>>);

}

#endif
