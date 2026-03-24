#ifndef HPP_GUARD_CTRLPP_ESTIMATION_MEKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_MEKF_H

// Multiplicative Extended Kalman Filter (MEKF) for attitude estimation.
//
// Two-track state: 7D nominal (quaternion + bias) with 6D tangent-space
// error-state covariance. The mandatory post-update covariance reset via
// frame-change Jacobian G is the key correctness concern.
//
// References:
//   - Markley, "Attitude Error Representations for Kalman Filtering" (2003)
//   - Sola, "Quaternion kinematics for the error-state Kalman filter" (2017)

#include "ctrlpp/so3.h"
#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"

#include "ctrlpp/detail/numerical_mekf_diff.h"

#include <Eigen/Geometry>

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace ctrlpp {

template <typename M, typename Scalar, std::size_t NB, std::size_t NY>
concept mekf_measurement_model =
    requires(const M &m, const Eigen::Quaternion<Scalar> &q,
             const Vector<Scalar, NB> &b)
    {
        { m(q, b) } -> std::convertible_to<Vector<Scalar, NY>>;
    };

template <typename M, typename Scalar, std::size_t NB, std::size_t NY>
concept differentiable_mekf_measurement =
    mekf_measurement_model<M, Scalar, NB, NY> &&
    requires(const M &m, const Eigen::Quaternion<Scalar> &q,
             const Vector<Scalar, NB> &b)
    {
        { m.jacobian(q, b) } -> std::convertible_to<Matrix<Scalar, NY, 3 + NB>>;
    };

template <typename Scalar, std::size_t NB, std::size_t NY>
struct mekf_config
{
    static constexpr std::size_t NE = 3 + NB;
    Matrix<Scalar, NE, NE> Q{Matrix<Scalar, NE, NE>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
    Vector<Scalar, NB> b0{Vector<Scalar, NB>::Zero()};
    Matrix<Scalar, NE, NE> P0{Matrix<Scalar, NE, NE>::Identity()};
    Scalar dt{Scalar{0.01}};
    Scalar numerical_eps{std::sqrt(std::numeric_limits<Scalar>::epsilon())};
};

template <typename Scalar, std::size_t NB, std::size_t NY, typename Measurement>
    requires mekf_measurement_model<Measurement, Scalar, NB, NY>
class mekf
{
    static constexpr std::size_t NE = 3 + NB;
    static constexpr int ne = static_cast<int>(NE);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr int nb = static_cast<int>(NB);

public:
    using observer_tag = struct mekf_tag;
    using state_vector_t = Vector<Scalar, 4 + NB>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NE, NE>;
    using meas_cov_t = Matrix<Scalar, NY, NY>;

    mekf(Measurement measurement, mekf_config<Scalar, NB, NY> config)
        : measurement_{std::move(measurement)}
        , q_{config.q0.normalized()}
        , b_{std::move(config.b0)}
        , P_{std::move(config.P0)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , eps_{config.numerical_eps}
        , dt_{config.dt}
        , innovation_{output_vector_t::Zero()}
        , state_cache_{}
    {
        update_state_cache();
    }

    void predict(const input_vector_t &omega)
    {
        predict_impl(omega, dt_);
    }

    void predict(const input_vector_t &omega, Scalar dt)
    {
        predict_impl(omega, dt);
    }

    void update(const output_vector_t &z)
    {
        // Predicted measurement
        auto z_pred = measurement_(q_, b_);
        innovation_ = (z - z_pred).eval();

        // Measurement Jacobian H (NY x NE)
        Matrix<Scalar, NY, NE> H;
        if constexpr(differentiable_mekf_measurement<Measurement, Scalar, NB, NY>)
        {
            H = measurement_.jacobian(q_, b_);
        }
        else
        {
            H = detail::numerical_mekf_jacobian<Scalar, NB, NY>(measurement_, q_, b_, eps_);
        }

        // Innovation covariance
        meas_cov_t S = (H * P_ * H.transpose() + R_).eval();

        // Kalman gain via transpose-solve
        Eigen::Matrix<Scalar, ny, ne> KT = S.transpose().colPivHouseholderQr().solve(H * P_);
        Eigen::Matrix<Scalar, ne, ny> K = KT.transpose().eval();

        // Error-state correction
        Vector<Scalar, NE> delta_xi = K * innovation_;

        // Multiplicative quaternion injection
        q_ = (q_ * so3::exp(delta_xi.template head<3>())).normalized();
        b_ += delta_xi.template tail<NB>();

        // Joseph-form covariance update
        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        P_ = (IKH * P_ * IKH.transpose() + K * R_ * K.transpose()).eval();

        // Mandatory covariance reset via frame-change Jacobian G
        cov_matrix_t G = cov_matrix_t::Identity();
        G.template block<3, 3>(0, 0) -= Scalar{0.5} * so3::skew(delta_xi.template head<3>());
        P_ = (G * P_ * G.transpose()).eval();

        // Symmetrize
        P_ = (Scalar{0.5} * (P_ + P_.transpose())).eval();

        update_state_cache();
    }

    [[nodiscard]] auto state() const -> const state_vector_t & { return state_cache_; }
    [[nodiscard]] auto covariance() const -> const cov_matrix_t & { return P_; }
    [[nodiscard]] auto innovation() const -> const output_vector_t & { return innovation_; }
    [[nodiscard]] auto attitude() const -> Eigen::Quaternion<Scalar> { return q_; }
    [[nodiscard]] auto bias() const -> const Vector<Scalar, NB> & { return b_; }

private:
    void predict_impl(const input_vector_t &omega, Scalar dt)
    {
        // Corrected angular velocity (remove estimated bias)
        Vector<Scalar, 3> omega_corr = omega - b_.template head<3>();

        // Nominal quaternion propagation
        q_ = (q_ * so3::exp(omega_corr * dt)).normalized();

        // Rotation matrix from corrected angular increment
        Eigen::Matrix<Scalar, 3, 3> C = so3::exp(omega_corr * dt).toRotationMatrix();

        // Error-state transition matrix F (NE x NE)
        cov_matrix_t F = cov_matrix_t::Identity();
        F.template block<3, 3>(0, 0) = C;
        F.template block<3, nb>(0, 3) = -Eigen::Matrix<Scalar, 3, nb>::Identity() * dt;

        // Propagate covariance
        P_ = (F * P_ * F.transpose() + Q_).eval();
        P_ = (Scalar{0.5} * (P_ + P_.transpose())).eval();

        update_state_cache();
    }

    void update_state_cache()
    {
        state_cache_(0) = q_.w();
        state_cache_(1) = q_.x();
        state_cache_(2) = q_.y();
        state_cache_(3) = q_.z();
        for(std::size_t i = 0; i < NB; ++i)
        {
            state_cache_(static_cast<Eigen::Index>(4 + i)) = b_(static_cast<Eigen::Index>(i));
        }
    }

    Measurement measurement_;
    Eigen::Quaternion<Scalar> q_;
    Vector<Scalar, NB> b_;
    cov_matrix_t P_;
    cov_matrix_t Q_;
    meas_cov_t R_;
    Scalar eps_;
    Scalar dt_;
    output_vector_t innovation_;
    state_vector_t state_cache_;
};

// CTAD deduction guide
template <typename Measurement, typename Scalar, std::size_t NB, std::size_t NY>
mekf(Measurement, mekf_config<Scalar, NB, NY>)
    -> mekf<Scalar, NB, NY, Measurement>;

namespace detail {

struct mekf_sa_measurement
{
    auto operator()(const Eigen::Quaternion<double> &,
                    const Vector<double, 3> &) const -> Vector<double, 3>
    {
        return Vector<double, 3>::Zero();
    }
};

}

static_assert(ObserverPolicy<mekf<double, 3, 3, detail::mekf_sa_measurement>>);
static_assert(CovarianceObserver<mekf<double, 3, 3, detail::mekf_sa_measurement>>);

}

#endif
