#ifndef HPP_GUARD_CTRLPP_UKF_H
#define HPP_GUARD_CTRLPP_UKF_H

#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"
#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/measurement_model.h"
#include "ctrlpp/sigma_points/sigma_point_strategy.h"
#include "ctrlpp/sigma_points/merwe_sigma_points.h"

#include <Eigen/Cholesky>

#include <cstddef>
#include <utility>

namespace ctrlpp {

enum class gain_decomposition { ldlt, qr };

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ukf_config {
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    gain_decomposition decomposition{gain_decomposition::ldlt};
};

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         typename Dynamics, typename Measurement,
         typename Strategy = merwe_sigma_points<Scalar, NX>>
requires dynamics_model<Dynamics, Scalar, NX, NU> &&
         measurement_model<Measurement, Scalar, NX, NY> &&
         sigma_point_strategy<Strategy, Scalar, NX>
class ukf {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr std::size_t num_sigma = Strategy::num_points;

public:
    using observer_tag = struct ukf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;

    ukf(Dynamics dynamics, Measurement measurement,
        ukf_config<Scalar, NX, NU, NY> config)
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , x_{std::move(config.x0)}
        , P_{std::move(config.P0)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , decomposition_{config.decomposition}
        , strategy_{}
        , innovation_{output_vector_t::Zero()}
    {
    }

    ukf(Dynamics dynamics, Measurement measurement,
        ukf_config<Scalar, NX, NU, NY> config,
        typename Strategy::options_t strategy_options)
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , x_{std::move(config.x0)}
        , P_{std::move(config.P0)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , decomposition_{config.decomposition}
        , strategy_{std::move(strategy_options)}
        , innovation_{output_vector_t::Zero()}
    {
    }

    ukf(Dynamics dynamics, Measurement measurement,
        ukf_config<Scalar, NX, NU, NY> config,
        Strategy strategy)
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , x_{std::move(config.x0)}
        , P_{std::move(config.P0)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , decomposition_{config.decomposition}
        , strategy_{std::move(strategy)}
        , innovation_{output_vector_t::Zero()}
    {
    }

    void predict(const input_vector_t& u)
    {
        auto sigma = strategy_.generate(x_, P_);

        // Propagate sigma points through dynamics
        std::array<state_vector_t, num_sigma> propagated;
        for (std::size_t i = 0; i < num_sigma; ++i) {
            propagated[i] = dynamics_(sigma.points[i], u);
        }

        // Weighted mean
        state_vector_t x_pred = state_vector_t::Zero();
        for (std::size_t i = 0; i < num_sigma; ++i) {
            x_pred += sigma.Wm[i] * propagated[i];
        }

        // Weighted covariance + Q
        cov_matrix_t P_pred = cov_matrix_t::Zero();
        for (std::size_t i = 0; i < num_sigma; ++i) {
            auto diff = (propagated[i] - x_pred).eval();
            P_pred += sigma.Wc[i] * diff * diff.transpose();
        }
        P_pred += Q_;

        // Symmetrize
        x_ = x_pred;
        P_ = (Scalar{0.5} * (P_pred + P_pred.transpose())).eval();
    }

    void update(const output_vector_t& z)
    {
        // Regenerate sigma points from predicted state
        auto sigma = strategy_.generate(x_, P_);

        // Transform sigma points through measurement model
        std::array<output_vector_t, num_sigma> z_sigma;
        for (std::size_t i = 0; i < num_sigma; ++i) {
            z_sigma[i] = measurement_(sigma.points[i]);
        }

        // Predicted measurement (weighted mean)
        output_vector_t z_pred = output_vector_t::Zero();
        for (std::size_t i = 0; i < num_sigma; ++i) {
            z_pred += sigma.Wm[i] * z_sigma[i];
        }

        // Innovation covariance S + R
        meas_cov_matrix_t S = meas_cov_matrix_t::Zero();
        for (std::size_t i = 0; i < num_sigma; ++i) {
            auto dz = (z_sigma[i] - z_pred).eval();
            S += sigma.Wc[i] * dz * dz.transpose();
        }
        S += R_;

        // Cross covariance Pxz
        Eigen::Matrix<Scalar, nx, ny> Pxz = Eigen::Matrix<Scalar, nx, ny>::Zero();
        for (std::size_t i = 0; i < num_sigma; ++i) {
            auto dx = (sigma.points[i] - x_).eval();
            auto dz = (z_sigma[i] - z_pred).eval();
            Pxz += sigma.Wc[i] * dx * dz.transpose();
        }

        // Kalman gain K = Pxz * S^-1
        Eigen::Matrix<Scalar, nx, ny> K;
        if (decomposition_ == gain_decomposition::ldlt) {
            // K = Pxz * S^-1  =>  S^T * K^T = Pxz^T  =>  solve for K^T
            Eigen::Matrix<Scalar, ny, nx> KT =
                S.ldlt().solve(Pxz.transpose());
            K = KT.transpose();
        } else {
            Eigen::Matrix<Scalar, ny, nx> KT =
                S.transpose().colPivHouseholderQr().solve(Pxz.transpose());
            K = KT.transpose();
        }

        // Innovation
        innovation_ = (z - z_pred).eval();

        // State update
        x_ = (x_ + K * innovation_).eval();

        // Covariance update: P -= K * S * K^T
        cov_matrix_t P_new = (P_ - K * S * K.transpose()).eval();

        // Symmetrize
        P_ = (Scalar{0.5} * (P_new + P_new.transpose())).eval();
    }

    [[nodiscard]] auto state() const -> const state_vector_t& { return x_; }
    [[nodiscard]] auto covariance() const -> const cov_matrix_t& { return P_; }
    [[nodiscard]] auto innovation() const -> const output_vector_t& { return innovation_; }

private:
    Dynamics dynamics_;
    Measurement measurement_;
    state_vector_t x_;
    cov_matrix_t P_;
    cov_matrix_t Q_;
    meas_cov_matrix_t R_;
    gain_decomposition decomposition_;
    Strategy strategy_;
    output_vector_t innovation_;
};

// CTAD deduction guide
template<typename Dynamics, typename Measurement, typename Scalar,
         std::size_t NX, std::size_t NU, std::size_t NY>
ukf(Dynamics, Measurement, ukf_config<Scalar, NX, NU, NY>)
    -> ukf<Scalar, NX, NU, NY, Dynamics, Measurement, merwe_sigma_points<Scalar, NX>>;

namespace detail {

struct ukf_sa_dynamics {
    auto operator()(const Vector<double, 2>&,
                    const Vector<double, 1>&) const -> Vector<double, 2>
    {
        return Vector<double, 2>::Zero();
    }
};

struct ukf_sa_measurement {
    auto operator()(const Vector<double, 2>&) const -> Vector<double, 1>
    {
        return Vector<double, 1>::Zero();
    }
};

}

static_assert(ObserverPolicy<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);

}

#endif
