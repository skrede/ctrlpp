#ifndef HPP_GUARD_CTRLPP_ESTIMATION_UKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_UKF_H

/// @brief Unscented Kalman Filter with swappable sigma point strategies.
///
/// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001

#include "ctrlpp/types.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <Eigen/Cholesky>

#include <cstddef>
#include <utility>

namespace ctrlpp
{

enum class gain_decomposition
{
    ldlt,
    qr
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ukf_config
{
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    gain_decomposition decomposition{gain_decomposition::ldlt};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy = merwe_sigma_points<Scalar, NX>>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && sigma_point_strategy<Strategy, Scalar, NX>
class ukf
{
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

    ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_x{std::move(config.x0)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_decomposition{config.decomposition}
        , m_strategy{}
        , m_innovation{output_vector_t::Zero()}
    {
    }

    ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, typename Strategy::options_t strategy_options)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_x{std::move(config.x0)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_decomposition{config.decomposition}
        , m_strategy{std::move(strategy_options)}
        , m_innovation{output_vector_t::Zero()}
    {
    }

    ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, Strategy strategy)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_x{std::move(config.x0)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_decomposition{config.decomposition}
        , m_strategy{std::move(strategy)}
        , m_innovation{output_vector_t::Zero()}
    {
    }

    void predict(const input_vector_t& u)
    {
        auto sigma = m_strategy.generate(m_x, m_P);

        // Propagate sigma points through dynamics
        std::array<state_vector_t, num_sigma> propagated;
        for(std::size_t i = 0; i < num_sigma; ++i)
            propagated[i] = m_dynamics(sigma.points[i], u);

        // Weighted mean
        state_vector_t x_pred = state_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            x_pred += sigma.Wm[i] * propagated[i];

        // Weighted covariance + Q
        cov_matrix_t P_pred = cov_matrix_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto diff = (propagated[i] - x_pred).eval();
            P_pred += sigma.Wc[i] * diff * diff.transpose();
        }
        P_pred += m_Q;

        // Symmetrize
        m_x = x_pred;
        m_P = (Scalar{0.5} * (P_pred + P_pred.transpose())).eval();
    }

    void update(const output_vector_t& z)
    {
        // Regenerate sigma points from predicted state
        auto sigma = m_strategy.generate(m_x, m_P);

        // Transform sigma points through measurement model
        std::array<output_vector_t, num_sigma> z_sigma;
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_sigma[i] = m_measurement(sigma.points[i]);

        // Predicted measurement (weighted mean)
        output_vector_t z_pred = output_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_pred += sigma.Wm[i] * z_sigma[i];

        // Innovation covariance S + R
        meas_cov_matrix_t S = meas_cov_matrix_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto dz = (z_sigma[i] - z_pred).eval();
            S += sigma.Wc[i] * dz * dz.transpose();
        }
        S += m_R;

        // Cross covariance Pxz
        Eigen::Matrix<Scalar, nx, ny> Pxz = Eigen::Matrix<Scalar, nx, ny>::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto dx = (sigma.points[i] - m_x).eval();
            auto dz = (z_sigma[i] - z_pred).eval();
            Pxz += sigma.Wc[i] * dx * dz.transpose();
        }

        // Kalman gain K = Pxz * S^-1
        Eigen::Matrix<Scalar, nx, ny> K;
        if(m_decomposition == gain_decomposition::ldlt)
        {
            // K = Pxz * S^-1  =>  S^T * K^T = Pxz^T  =>  solve for K^T
            Eigen::Matrix<Scalar, ny, nx> KT = S.ldlt().solve(Pxz.transpose());
            K = KT.transpose();
        }
        else
        {
            Eigen::Matrix<Scalar, ny, nx> KT = S.transpose().colPivHouseholderQr().solve(Pxz.transpose());
            K = KT.transpose();
        }

        // Innovation
        m_innovation = (z - z_pred).eval();

        // State update
        m_x = (m_x + K * m_innovation).eval();

        // Covariance update: P -= K * S * K^T
        cov_matrix_t P_new = (m_P - K * S * K.transpose()).eval();

        // Symmetrize
        m_P = (Scalar{0.5} * (P_new + P_new.transpose())).eval();
    }

    const state_vector_t& state() const { return m_x; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

private:
    Dynamics m_dynamics;
    Measurement m_measurement;
    state_vector_t m_x;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    meas_cov_matrix_t m_R;
    gain_decomposition m_decomposition;
    Strategy m_strategy;
    output_vector_t m_innovation;
};

// CTAD deduction guide
template <typename Dynamics, typename Measurement, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
ukf(Dynamics, Measurement, ukf_config<Scalar, NX, NU, NY>) -> ukf<Scalar, NX, NU, NY, Dynamics, Measurement, merwe_sigma_points<Scalar, NX>>;

namespace detail
{

struct ukf_sa_dynamics
{
    auto operator()(const Vector<double, 2>&, const Vector<double, 1>&) const -> Vector<double, 2> { return Vector<double, 2>::Zero(); }
};

struct ukf_sa_measurement
{
    auto operator()(const Vector<double, 2>&) const -> Vector<double, 1> { return Vector<double, 1>::Zero(); }
};

} // namespace detail

static_assert(ObserverPolicy<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);

} // namespace ctrlpp

#endif
