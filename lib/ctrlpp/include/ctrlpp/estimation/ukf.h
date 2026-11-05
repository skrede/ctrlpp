#ifndef HPP_GUARD_CTRLPP_ESTIMATION_UKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_UKF_H

// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001

#include "ctrlpp/estimation/ukf_types.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include "ctrlpp/estimation/detail/filter_configuration.h"

#include <Eigen/Cholesky>

#include <array>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

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

    struct validated_tag
    {
    };

public:
    using observer_tag = struct ukf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;

    static ctrlpp::expected<ukf, filter_error> create(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, Strategy strategy = Strategy{})
    {
        if(const auto valid = detail::validate_filter_configuration(config.Q, config.R, config.x0, config.P0); !valid)
            return ctrlpp::unexpected(valid.error());
        return ukf{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config), std::move(strategy)};
    }

    static ctrlpp::expected<ukf, filter_error> create(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, typename Strategy::options_t strategy_options)
    {
        auto strategy = Strategy::try_create(std::move(strategy_options));
        if(!strategy)
            return ctrlpp::unexpected(strategy.error());
        return create(std::move(dynamics), std::move(measurement), std::move(config), std::move(*strategy));
    }

    // Not fallible: the plant already took the command. A prediction that poisons the estimate or is not
    // positive definite is latched in health() instead.
    void predict(const input_vector_t& u)
    {
        auto sigma = m_strategy.generate(m_x, m_P);
        note_health(sigma);
        auto propagated = propagate_sigma_points(sigma.points, u);
        auto x_pred = compute_predicted_mean(sigma.Wm, propagated);
        m_P = compute_predicted_covariance(sigma.Wc, propagated, x_pred);
        m_x = x_pred;
        note_prediction();
    }

    // A refused step leaves every member bitwise unchanged, the health included, so the caller may retry
    // with the next sample.
    ctrlpp::expected<void, ukf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_x = m_x;
        auto const previous_P = m_P;
        auto const previous_innovation = m_innovation;
        auto const previous_nis = m_nis;
        auto const previous_health = m_health;
        auto sigma = m_strategy.generate(m_x, m_P);
        note_health(sigma);
        auto z_sigma = compute_measurement_sigma_points(sigma.points);
        auto z_pred = compute_predicted_measurement(sigma.Wm, z_sigma);
        auto [S, Pxz] = compute_innovation_and_cross_covariance(sigma, z_sigma, z_pred);
        auto K = compute_kalman_gain(Pxz, S);

        m_innovation = (z - z_pred).eval();
        m_nis = (m_innovation.transpose() * S.colPivHouseholderQr().solve(m_innovation))(0, 0);
        apply_correction_and_update_covariance(K, S, reduction_base(sigma));
        if(!m_x.allFinite() || !m_P.allFinite() || !m_innovation.allFinite() || !std::isfinite(m_nis))
            return roll_back(previous_x, previous_P, previous_innovation, previous_nis, previous_health);
        note_definiteness();
        return {};
    }

    const state_vector_t& state() const { return m_x; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

    Scalar nis() const { return m_nis; }

    ukf_health health() const { return m_health; }

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
    Scalar m_nis;
    ukf_health m_health;

    ukf(validated_tag, Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config, Strategy strategy)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_x{std::move(config.x0)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_decomposition{config.decomposition}
        , m_strategy{std::move(strategy)}
        , m_innovation{output_vector_t::Zero()}
        , m_nis{0}
        , m_health{ukf_health::ok}
    {
    }

    ctrlpp::expected<void, ukf_update_error> roll_back(const state_vector_t& x, const cov_matrix_t& P, const output_vector_t& innovation, Scalar nis, ukf_health health)
    {
        m_x = x;
        m_P = P;
        m_innovation = innovation;
        m_nis = nis;
        m_health = health;
        return ctrlpp::unexpected(ukf_update_error::non_finite_result);
    }

    // Classifies before anything is mutated, the carried estimate ahead of the measurement.
    ctrlpp::expected<void, ukf_update_error> check_step(const output_vector_t& z) const
    {
        if(!m_x.allFinite())
            return ctrlpp::unexpected(ukf_update_error::non_finite_state);
        if(!m_P.allFinite())
            return ctrlpp::unexpected(ukf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(ukf_update_error::non_finite_measurement);
        return {};
    }

    ukf_update_error latch_health(ukf_update_error fault)
    {
        if(fault != ukf_update_error::non_finite_measurement)
            m_health = ukf_health::non_finite_estimate;
        return fault;
    }

    template <typename Sigma>
    void note_health(const Sigma& sigma)
    {
        if(sigma.spd_repaired && m_health == ukf_health::ok)
            m_health = ukf_health::covariance_repaired;
    }

    void note_prediction()
    {
        if(!m_x.allFinite() || !m_P.allFinite())
            m_health = ukf_health::non_finite_estimate;
        else
            note_definiteness();
    }

    // Judged by the unpivoted Cholesky factorization the Merwe and Julier generations attempt first on this
    // same matrix, so the status never disagrees with the repair the next generation makes.
    void note_definiteness()
    {
        if(m_health == ukf_health::ok && Eigen::LLT<cov_matrix_t>{m_P}.info() != Eigen::Success)
            m_health = ukf_health::covariance_repaired;
    }

    template <typename Sigma>
    cov_matrix_t reduction_base(const Sigma& sigma) const
    {
        if(sigma.spd_repaired)
            return (sigma.factor * sigma.factor.transpose()).eval();
        return m_P;
    }

    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Sec. 3.1
    std::array<state_vector_t, num_sigma> propagate_sigma_points(const std::array<state_vector_t, num_sigma>& points, const input_vector_t& u) const
    {
        std::array<state_vector_t, num_sigma> propagated;
        for(std::size_t i = 0; i < num_sigma; ++i)
            propagated[i] = m_dynamics(points[i], u);
        return propagated;
    }

    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 17
    state_vector_t compute_predicted_mean(const std::array<Scalar, num_sigma>& Wm, const std::array<state_vector_t, num_sigma>& propagated) const
    {
        state_vector_t x_pred = state_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            x_pred += Wm[i] * propagated[i];
        return x_pred;
    }

    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 18
    cov_matrix_t compute_predicted_covariance(const std::array<Scalar, num_sigma>& Wc, const std::array<state_vector_t, num_sigma>& propagated, const state_vector_t& x_pred) const
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

    std::array<output_vector_t, num_sigma> compute_measurement_sigma_points(const std::array<state_vector_t, num_sigma>& points) const
    {
        std::array<output_vector_t, num_sigma> z_sigma;
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_sigma[i] = m_measurement(points[i]);
        return z_sigma;
    }

    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 22
    output_vector_t compute_predicted_measurement(const std::array<Scalar, num_sigma>& Wm, const std::array<output_vector_t, num_sigma>& z_sigma) const
    {
        output_vector_t z_pred = output_vector_t::Zero();
        for(std::size_t i = 0; i < num_sigma; ++i)
            z_pred += Wm[i] * z_sigma[i];
        return z_pred;
    }

    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 23-24
    std::pair<meas_cov_matrix_t, Eigen::Matrix<Scalar, nx, ny>> compute_innovation_and_cross_covariance(
            const auto& sigma, const std::array<output_vector_t, num_sigma>& z_sigma, const output_vector_t& z_pred) const
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

    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 25
    Eigen::Matrix<Scalar, nx, ny> compute_kalman_gain(const Eigen::Matrix<Scalar, nx, ny>& Pxz, const meas_cov_matrix_t& S) const
    {
        Eigen::Matrix<Scalar, ny, nx> KT;
        if(m_decomposition == gain_decomposition::ldlt)
            KT = S.ldlt().solve(Pxz.transpose());
        else if constexpr(ny == 1)
            KT = detail::one_output_gain_transpose(S(0, 0), Pxz.transpose());
        else
            KT = S.transpose().colPivHouseholderQr().solve(Pxz.transpose());
        return KT.transpose().eval();
    }

    // The reduction P - K S K^T equals P - K Pxz^T, so no K R K^T term is added. It starts from the
    // covariance the sigma set represents, the repaired one after a repair, since K and S describe that one.
    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 26-27
    void apply_correction_and_update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K, const meas_cov_matrix_t& S, const cov_matrix_t& base)
    {
        m_x = (m_x + K * m_innovation).eval();
        m_P = detail::symmetrize((base - K * S * K.transpose()).eval());
    }
};

template <typename Dynamics, typename Measurement, ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
ukf(Dynamics, Measurement, ukf_config<Scalar, NX, NU, NY>) -> ukf<Scalar, NX, NU, NY, Dynamics, Measurement, merwe_sigma_points<Scalar, NX>>;

namespace detail
{

struct ukf_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>&, const Vector<double, 1>&) const { return Vector<double, 2>::Zero(); }
};

struct ukf_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>&) const { return Vector<double, 1>::Zero(); }
};

}

static_assert(ObserverPolicy<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, detail::ukf_sa_dynamics, detail::ukf_sa_measurement>>);

}

#endif
