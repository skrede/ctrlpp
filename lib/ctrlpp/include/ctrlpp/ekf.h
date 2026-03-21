#ifndef HPP_GUARD_CTRLPP_EKF_H
#define HPP_GUARD_CTRLPP_EKF_H

#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"
#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/measurement_model.h"
#include "ctrlpp/mpc/differentiable_dynamics.h"
#include "ctrlpp/mpc/differentiable_measurement.h"
#include "ctrlpp/detail/numerical_diff.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ekf_config {
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar numerical_eps{std::sqrt(std::numeric_limits<Scalar>::epsilon())};
};

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         typename Dynamics, typename Measurement>
requires dynamics_model<Dynamics, Scalar, NX, NU> &&
         measurement_model<Measurement, Scalar, NX, NY>
class ekf {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct ekf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;

    ekf(Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config)
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , x_{std::move(config.x0)}
        , P_{std::move(config.P0)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , eps_{config.numerical_eps}
        , innovation_{output_vector_t::Zero()}
    {
    }

    void predict(const input_vector_t& u)
    {
        auto x_prev = x_;

        x_ = dynamics_(x_prev, u);

        Matrix<Scalar, NX, NX> F;
        if constexpr (differentiable_dynamics<Dynamics, Scalar, NX, NU>) {
            F = dynamics_.jacobian_x(x_prev, u);
        } else {
            F = detail::numerical_jacobian_x<Scalar, NX, NU>(dynamics_, x_prev, u, eps_);
        }

        P_ = (F * P_ * F.transpose() + Q_).eval();
        P_ = (Scalar{0.5} * (P_ + P_.transpose())).eval();
    }

    void update(const output_vector_t& z)
    {
        auto z_pred = measurement_(x_);

        innovation_ = (z - z_pred).eval();

        Matrix<Scalar, NY, NX> H;
        if constexpr (differentiable_measurement<Measurement, Scalar, NX, NY>) {
            H = measurement_.jacobian(x_);
        } else {
            H = detail::numerical_jacobian_h<Scalar, NX, NY>(measurement_, x_, eps_);
        }

        meas_cov_matrix_t S = (H * P_ * H.transpose() + R_).eval();

        Eigen::Matrix<Scalar, ny, nx> KT_solved =
            S.transpose().colPivHouseholderQr().solve(H * P_);
        Eigen::Matrix<Scalar, nx, ny> K = KT_solved.transpose().eval();

        x_ = (x_ + K * innovation_).eval();

        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        P_ = (IKH * P_ * IKH.transpose() + K * R_ * K.transpose()).eval();
        P_ = (Scalar{0.5} * (P_ + P_.transpose())).eval();

        nees_ = (innovation_.transpose() * S.colPivHouseholderQr().solve(innovation_))(0, 0);
    }

    [[nodiscard]] auto state() const -> const state_vector_t& { return x_; }
    [[nodiscard]] auto covariance() const -> const cov_matrix_t& { return P_; }
    [[nodiscard]] auto innovation() const -> const output_vector_t& { return innovation_; }
    [[nodiscard]] auto nees() const -> Scalar { return nees_; }

private:
    Dynamics dynamics_;
    Measurement measurement_;
    state_vector_t x_;
    cov_matrix_t P_;
    cov_matrix_t Q_;
    meas_cov_matrix_t R_;
    Scalar eps_;
    output_vector_t innovation_;
    Scalar nees_{0};
};

// CTAD deduction guide
template<typename Dynamics, typename Measurement, typename Scalar,
         std::size_t NX, std::size_t NU, std::size_t NY>
ekf(Dynamics, Measurement, ekf_config<Scalar, NX, NU, NY>)
    -> ekf<Scalar, NX, NU, NY, Dynamics, Measurement>;

namespace detail {

struct ekf_sa_dynamics {
    auto operator()(const Vector<double, 2>&,
                    const Vector<double, 1>&) const -> Vector<double, 2>
    {
        return Vector<double, 2>::Zero();
    }
};

struct ekf_sa_measurement {
    auto operator()(const Vector<double, 2>&) const -> Vector<double, 1>
    {
        return Vector<double, 1>::Zero();
    }
};

}

static_assert(ObserverPolicy<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);
static_assert(CovarianceObserver<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);

}

#endif
