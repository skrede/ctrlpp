#ifndef HPP_GUARD_CTRLPP_ESTIMATION_UKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_UKF_H

/// @brief Unscented Kalman Filter with swappable sigma point strategies.
///
/// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001

#include "ctrlpp/types.h"
#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <Eigen/Cholesky>

#include <array>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

enum class gain_decomposition
{
    ldlt,
    qr
};

/// @brief Filter-health status the caller can inspect after any step.
///
/// `ok` means every covariance factored so far was positive definite;
/// `covariance_repaired` latches once a non-positive-definite covariance had
/// to be repaired to the nearest symmetric positive definite matrix, signaling
/// that the estimate has entered a numerically degraded regime.
enum class ukf_health
{
    ok,
    covariance_repaired
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ukf_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    gain_decomposition decomposition{gain_decomposition::ldlt};
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement, typename Strategy = merwe_sigma_points<Scalar, NX>>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && sigma_point_strategy<Strategy, Scalar, NX>
class ukf
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

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

    /// @brief Fallible factory for construction from a sigma-point strategy's
    /// options aggregate.
    ///
    /// This is the path that builds the strategy inside the filter, so it is
    /// where an out-of-domain parameter set would otherwise pass unreported:
    /// the strategy's own rejection is forwarded verbatim rather than swallowed.
    /// The overload taking an already-built strategy needs no such forwarding,
    /// since the strategy was validated where it was constructed.
    ///
    /// Requires the strategy to expose a `try_create` returning
    /// `ctrlpp::expected<Strategy, filter_error>`; a strategy without one has
    /// no failure to forward and is constructed directly instead.
    [[nodiscard]] static auto try_create(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, typename Strategy::options_t strategy_options)
        -> ctrlpp::expected<ukf, filter_error>
    {
        auto strategy = Strategy::try_create(std::move(strategy_options));
        if(!strategy)
            return ctrlpp::unexpected(strategy.error());
        return ukf{std::move(dynamics), std::move(measurement), std::move(config), std::move(*strategy)};
    }

    /// @brief Construct from a sigma-point strategy's options aggregate.
    ///
    /// The strategy is built inside the filter, so this constructor can only
    /// report an out-of-domain parameter set the way the strategy itself does.
    /// With a strategy whose options constructor is the exception-gated wrapper
    /// over its `try_create`, that means this overload throws in the exceptions
    /// tree and does not compile on an exception-free build, where
    /// `ukf::try_create` is the construction path to use.
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
        note_health(sigma);
        auto propagated = propagate_sigma_points(sigma.points, u);
        auto x_pred = compute_predicted_mean(sigma.Wm, propagated);
        m_P = compute_predicted_covariance(sigma.Wc, propagated, x_pred);
        m_x = x_pred;
    }

    void update(const output_vector_t& z)
    {
        auto sigma = m_strategy.generate(m_x, m_P);
        note_health(sigma);
        auto z_sigma = compute_measurement_sigma_points(sigma.points);
        auto z_pred = compute_predicted_measurement(sigma.Wm, z_sigma);
        auto [S, Pxz] = compute_innovation_and_cross_covariance(sigma, z_sigma, z_pred);
        auto K = compute_kalman_gain(Pxz, S);

        m_innovation = (z - z_pred).eval();
        m_nis = (m_innovation.transpose() * S.colPivHouseholderQr().solve(m_innovation))(0, 0);
        apply_correction_and_update_covariance(K, S);
    }

    const state_vector_t& state() const { return m_x; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

    /// @brief Normalized Innovation Squared: innovation^T S^{-1} innovation
    /// (chi-square distributed with dof = NY under a consistent filter).
    Scalar nis() const { return m_nis; }

    /// @brief Report whether the filter has had to repair a non-positive-definite
    /// covariance to the nearest symmetric positive definite matrix.
    ukf_health health() const { return m_health; }

private:
    /// @brief Latch a degraded health status when a sigma-point generation had
    /// to repair a non-positive-definite covariance.
    template <typename Sigma>
    void note_health(const Sigma& sigma)
    {
        if(sigma.spd_repaired)
            m_health = ukf_health::covariance_repaired;
    }

    /// @brief Propagate sigma points through dynamics model.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Sec. 3.1
    auto propagate_sigma_points(const std::array<state_vector_t, num_sigma>& points, const input_vector_t& u) const -> std::array<state_vector_t, num_sigma>
    {
        std::array<state_vector_t, num_sigma> propagated;
        for(std::size_t i = 0; i < num_sigma; ++i)
            propagated[i] = m_dynamics(points[i], u);
        return propagated;
    }

    /// @brief Compute weighted mean of propagated sigma points.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 17
    auto compute_predicted_mean(const std::array<Scalar, num_sigma>& Wm, const std::array<state_vector_t, num_sigma>& propagated) const -> state_vector_t
    {
        state_vector_t x_pred = state_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            x_pred += Wm[i] * propagated[i];
        return x_pred;
    }

    /// @brief Compute weighted covariance of propagated sigma points + process noise Q.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 18
    auto compute_predicted_covariance(const std::array<Scalar, num_sigma>& Wc, const std::array<state_vector_t, num_sigma>& propagated, const state_vector_t& x_pred) const -> cov_matrix_t
    {
        cov_matrix_t P_pred = cov_matrix_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto diff = (propagated[i] - x_pred).eval();
            P_pred += Wc[i] * diff * diff.transpose();
        }
        P_pred += m_Q;
        return detail::symmetrize(P_pred);
    }

    /// @brief Transform sigma points through measurement model.
    auto compute_measurement_sigma_points(const std::array<state_vector_t, num_sigma>& points) const -> std::array<output_vector_t, num_sigma>
    {
        std::array<output_vector_t, num_sigma> z_sigma;
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_sigma[i] = m_measurement(points[i]);
        return z_sigma;
    }

    /// @brief Compute weighted mean of measurement sigma points.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 22
    auto compute_predicted_measurement(const std::array<Scalar, num_sigma>& Wm, const std::array<output_vector_t, num_sigma>& z_sigma) const -> output_vector_t
    {
        output_vector_t z_pred = output_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_pred += Wm[i] * z_sigma[i];
        return z_pred;
    }

    /// @brief Compute innovation covariance S and cross-covariance Pxz.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 23-24
    auto compute_innovation_and_cross_covariance(const auto& sigma, const std::array<output_vector_t, num_sigma>& z_sigma, const output_vector_t& z_pred) const
        -> std::pair<meas_cov_matrix_t, Eigen::Matrix<Scalar, nx, ny>>
    {
        meas_cov_matrix_t S = meas_cov_matrix_t::Zero();
        Eigen::Matrix<Scalar, nx, ny> Pxz = Eigen::Matrix<Scalar, nx, ny>::Zero();

        for(std::size_t i = 0; i < num_sigma; ++i)
        {
            auto dz = (z_sigma[i] - z_pred).eval();
            S += sigma.Wc[i] * dz * dz.transpose();

            auto dx = (sigma.points[i] - m_x).eval();
            Pxz += sigma.Wc[i] * dx * dz.transpose();
        }
        S += m_R;

        return {S, Pxz};
    }

    /// @brief Compute Kalman gain: K = Pxz * S^{-1} via LDLT or QR decomposition.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 25
    auto compute_kalman_gain(const Eigen::Matrix<Scalar, nx, ny>& Pxz, const meas_cov_matrix_t& S) const -> Eigen::Matrix<Scalar, nx, ny>
    {
        Eigen::Matrix<Scalar, ny, nx> KT;
        if(m_decomposition == gain_decomposition::ldlt)
            KT = S.ldlt().solve(Pxz.transpose());
        else
            KT = S.transpose().colPivHouseholderQr().solve(Pxz.transpose());
        return KT.transpose().eval();
    }

    /// @brief Apply the state correction and reduce the covariance.
    ///
    /// The posterior covariance is the minimum mean-square-error reduction
    /// P = P - K*S*K^T, with the innovation covariance S = Pzz + R and gain
    /// K = Pxz*S^{-1}. Since K*S*K^T = K*Pxz^T, this is exactly the amount of
    /// uncertainty the measurement removes; no extra term is added.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 26-27
    void apply_correction_and_update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K, const meas_cov_matrix_t& S)
    {
        m_x = (m_x + K * m_innovation).eval();
        m_P = detail::symmetrize((m_P - K * S * K.transpose()).eval());
    }

    Dynamics m_dynamics;
    Measurement m_measurement;
    state_vector_t m_x;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    meas_cov_matrix_t m_R;
    gain_decomposition m_decomposition;
    Strategy m_strategy;
    output_vector_t m_innovation;
    Scalar m_nis{0};
    ukf_health m_health{ukf_health::ok};
};

// CTAD deduction guide
template <typename Dynamics, typename Measurement, ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
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

}

static_assert(ObserverPolicy<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);

}

#endif
