#ifndef HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H

/// @brief Bootstrap SIR particle filter with ESS-adaptive resampling and roughening.

#include "ctrlpp/types.h"
#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/estimation/resampling/resampling_strategy.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstddef>
#include <numbers>
#include <utility>
#include <algorithm>
#include <type_traits>

namespace ctrlpp
{

enum class extraction_method
{
    weighted_mean,
    map
};

enum class weight_representation
{
    log,
    linear
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct pf_config
{
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar ess_threshold{Scalar{-1}};
    Scalar roughening_scale{Scalar{0.2}};
    extraction_method extraction{extraction_method::weighted_mean};
    weight_representation weights{weight_representation::log};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement, typename Resampler = systematic_resampling, typename Rng = std::mt19937_64>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && std::uniform_random_bit_generator<Rng> && resampling_strategy<Resampler, Rng, NP>
class particle_filter
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct pf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;

    particle_filter(Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{})
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_resampler{}
        , m_rng{std::move(rng)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_R_inv{m_R.colPivHouseholderQr().inverse()}
        , m_log_det_2piR{compute_log_det_2piR()}
        , m_ess_threshold{config.ess_threshold < Scalar{0} ? static_cast<Scalar>(NP) / Scalar{2} : config.ess_threshold}
        , m_roughening_scale{config.roughening_scale}
        , m_extraction{config.extraction}
        , m_weight_mode{config.weights}
    {
        initialize_particles(config.x0, config.P0);
    }

    particle_filter(Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Resampler resampler, Rng rng = Rng{})
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_resampler{std::move(resampler)}
        , m_rng{std::move(rng)}
        , m_Q{std::move(config.Q)}
        , m_R{std::move(config.R)}
        , m_R_inv{m_R.colPivHouseholderQr().inverse()}
        , m_log_det_2piR{compute_log_det_2piR()}
        , m_ess_threshold{config.ess_threshold < Scalar{0} ? static_cast<Scalar>(NP) / Scalar{2} : config.ess_threshold}
        , m_roughening_scale{config.roughening_scale}
        , m_extraction{config.extraction}
        , m_weight_mode{config.weights}
    {
        initialize_particles(config.x0, config.P0);
    }

    /// @brief Propagate all particles through dynamics with process noise.
    ///
    /// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Eq. 6
    void predict(const input_vector_t& u)
    {
        propagate_particles(u);
    }

    /// @brief Update particle weights with measurement and resample if ESS is low.
    ///
    /// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Alg. 1
    void update(const output_vector_t& z)
    {
        if(m_weight_mode == weight_representation::log)
            update_log(z);
        else
            update_linear(z);
    }

    state_vector_t state() const
    {
        if(m_extraction == extraction_method::map)
            return map_estimate();
        return weighted_mean();
    }

    [[nodiscard]] auto weighted_mean() const -> state_vector_t
    {
        state_vector_t mean = state_vector_t::Zero();
        if(m_weight_mode == weight_representation::log)
        {
            auto lin = log_to_linear();
            for(std::size_t i = 0; i < NP; ++i)
                mean += lin[i] * m_particles[i];
        }
        else
        {
            for(std::size_t i = 0; i < NP; ++i)
                mean += m_linear_weights[i] * m_particles[i];
        }
        return mean;
    }

    [[nodiscard]] auto map_estimate() const -> state_vector_t
    {
        std::size_t best = 0;
        if(m_weight_mode == weight_representation::log)
        {
            for(std::size_t i = 1; i < NP; ++i)
            {
                if(m_log_weights[i] > m_log_weights[best])
                    best = i;
            }
        }
        else
        {
            for(std::size_t i = 1; i < NP; ++i)
            {
                if(m_linear_weights[i] > m_linear_weights[best])
                    best = i;
            }
        }
        return m_particles[best];
    }

    [[nodiscard]] auto particles() const -> const std::array<state_vector_t, NP>& { return m_particles; }

private:
    Dynamics m_dynamics;
    Measurement m_measurement;
    Resampler m_resampler;
    Rng m_rng;

    Matrix<Scalar, NX, NX> m_Q;
    Matrix<Scalar, NY, NY> m_R;
    Matrix<Scalar, NY, NY> m_R_inv;
    Scalar m_log_det_2piR;
    Scalar m_ess_threshold;
    Scalar m_roughening_scale;
    extraction_method m_extraction;
    weight_representation m_weight_mode;

    std::array<state_vector_t, NP> m_particles;
    std::array<Scalar, NP> m_log_weights{};
    std::array<Scalar, NP> m_linear_weights{};

    // Cholesky factor of Q for process noise sampling
    Eigen::Matrix<Scalar, nx, nx> m_Q_L{};

    Scalar compute_log_det_2piR() const
    {
        Scalar log_det = m_R.colPivHouseholderQr().logAbsDeterminant();
        return static_cast<Scalar>(NY) * std::log(Scalar{2} * std::numbers::pi_v<Scalar>) + log_det;
    }

    void initialize_particles(const state_vector_t& x0, const Matrix<Scalar, NX, NX>& P0)
    {
        Eigen::LLT<Eigen::Matrix<Scalar, nx, nx>> llt_P0(P0);
        Eigen::Matrix<Scalar, nx, nx> sqrt_P0 = llt_P0.matrixL();

        Eigen::LLT<Eigen::Matrix<Scalar, nx, nx>> llt_Q(m_Q);
        m_Q_L = llt_Q.matrixL();

        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});

        for(std::size_t i = 0; i < NP; ++i)
        {
            state_vector_t noise;
            for(std::size_t d = 0; d < NX; ++d)
                noise(static_cast<int>(d)) = normal(m_rng);
            m_particles[i] = x0 + sqrt_P0 * noise;
        }

        Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
        Scalar lin_uniform = Scalar{1} / static_cast<Scalar>(NP);
        m_log_weights.fill(log_uniform);
        m_linear_weights.fill(lin_uniform);
    }

    state_vector_t sample_process_noise()
    {
        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});
        state_vector_t noise;
        for(std::size_t d = 0; d < NX; ++d)
            noise(static_cast<int>(d)) = normal(m_rng);
        return m_Q_L * noise;
    }

    /// @brief Propagate all particles through dynamics model with additive process noise.
    ///
    /// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Eq. 6
    void propagate_particles(const input_vector_t& u)
    {
        for(std::size_t i = 0; i < NP; ++i)
        {
            m_particles[i] = m_dynamics(m_particles[i], u);
            m_particles[i] += sample_process_noise();
        }
    }

    /// @brief Compute Gaussian log-likelihood of innovation.
    ///
    /// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 63
    Scalar log_likelihood(const output_vector_t& z, const output_vector_t& z_pred) const
    {
        output_vector_t innov = z - z_pred;
        Scalar mahal = (innov.transpose() * m_R_inv * innov)(0, 0);
        return Scalar{-0.5} * mahal - Scalar{0.5} * m_log_det_2piR;
    }

    /// @brief Compute log-weights for all particles given measurement.
    void compute_log_weights(const output_vector_t& z)
    {
        for(std::size_t i = 0; i < NP; ++i)
        {
            output_vector_t z_pred = m_measurement(m_particles[i]);
            m_log_weights[i] += log_likelihood(z, z_pred);
        }
    }

    /// @brief Normalize log-weights via log-sum-exp trick.
    ///
    /// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. III-A
    void normalize_log_weights()
    {
        Scalar max_log_w = *std::max_element(m_log_weights.begin(), m_log_weights.end());
        Scalar sum_exp = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
            sum_exp += std::exp(m_log_weights[i] - max_log_w);
        Scalar log_sum = max_log_w + std::log(sum_exp);
        for(auto& lw : m_log_weights)
            lw -= log_sum;
    }

    /// @brief Compute linear weights for all particles given measurement.
    void compute_linear_weights(const output_vector_t& z)
    {
        for(std::size_t i = 0; i < NP; ++i)
        {
            output_vector_t z_pred = m_measurement(m_particles[i]);
            Scalar ll = log_likelihood(z, z_pred);
            m_linear_weights[i] *= std::exp(ll);
        }
    }

    /// @brief Normalize linear weights.
    void normalize_linear_weights()
    {
        Scalar sum = Scalar{0};
        for(auto w : m_linear_weights)
            sum += w;
        if(sum > Scalar{0})
            for(auto& w : m_linear_weights)
                w /= sum;
    }

    /// @brief Compute Effective Sample Size from current weights.
    ///
    /// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 51
    Scalar compute_ess_from_log() const
    {
        Scalar sum_w2 = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
        {
            Scalar w = std::exp(m_log_weights[i]);
            sum_w2 += w * w;
        }
        return Scalar{1} / sum_w2;
    }

    /// @brief Compute Effective Sample Size from linear weights.
    Scalar compute_ess_from_linear() const
    {
        Scalar sum_w2 = Scalar{0};
        for(auto w : m_linear_weights)
            sum_w2 += w * w;
        return Scalar{1} / sum_w2;
    }

    /// @brief Resample particles and reset to uniform weights (log mode).
    void resample_log()
    {
        auto lin = log_to_linear();
        std::array<std::size_t, NP> indices{};
        m_resampler.resample(lin, indices, m_rng);
        reindex_particles(indices);
        apply_roughening();

        Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
        m_log_weights.fill(log_uniform);
    }

    /// @brief Resample particles and reset to uniform weights (linear mode).
    void resample_linear()
    {
        std::array<std::size_t, NP> indices{};
        m_resampler.resample(m_linear_weights, indices, m_rng);
        reindex_particles(indices);
        apply_roughening();

        Scalar uniform = Scalar{1} / static_cast<Scalar>(NP);
        m_linear_weights.fill(uniform);
    }

    /// @brief Log-weight update: compute weights, normalize, resample if ESS low.
    void update_log(const output_vector_t& z)
    {
        compute_log_weights(z);
        normalize_log_weights();
        if(compute_ess_from_log() < m_ess_threshold)
            resample_log();
    }

    /// @brief Linear-weight update: compute weights, normalize, resample if ESS low.
    void update_linear(const output_vector_t& z)
    {
        compute_linear_weights(z);
        normalize_linear_weights();
        if(compute_ess_from_linear() < m_ess_threshold)
            resample_linear();
    }

    std::array<Scalar, NP> log_to_linear() const
    {
        std::array<Scalar, NP> lin;
        for(std::size_t i = 0; i < NP; ++i)
            lin[i] = std::exp(m_log_weights[i]);
        return lin;
    }

    void reindex_particles(const std::array<std::size_t, NP>& indices)
    {
        std::array<state_vector_t, NP> resampled;
        for(std::size_t i = 0; i < NP; ++i)
            resampled[i] = m_particles[indices[i]];
        m_particles = resampled;
    }

    /// @brief Roughening: add jitter to resampled particles to prevent degeneracy.
    ///
    /// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Sec. 5
    void apply_roughening()
    {
        if(m_roughening_scale <= Scalar{0})
            return;

        Scalar np_factor = std::pow(static_cast<Scalar>(NP), Scalar{-1} / static_cast<Scalar>(NX));
        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});

        for(std::size_t d = 0; d < NX; ++d)
        {
            Scalar sigma = compute_roughening_sigma(d, np_factor);
            if(sigma > Scalar{0})
                apply_particle_jitter(d, sigma, normal);
        }
    }

    Scalar compute_roughening_sigma(std::size_t d, Scalar np_factor) const
    {
        int di = static_cast<int>(d);
        Scalar min_val = m_particles[0](di);
        Scalar max_val = min_val;
        for(std::size_t i = 1; i < NP; ++i)
        {
            Scalar v = m_particles[i](di);
            min_val = std::min(min_val, v);
            max_val = std::max(max_val, v);
        }
        Scalar sigma = m_roughening_scale * (max_val - min_val) * np_factor;

        // Fallback: if all particles collapsed, use process noise scale
        if(sigma <= Scalar{0})
            sigma = m_roughening_scale * std::sqrt(m_Q(di, di));
        return sigma;
    }

    void apply_particle_jitter(std::size_t d, Scalar sigma, std::normal_distribution<Scalar>& normal)
    {
        int di = static_cast<int>(d);
        for(std::size_t i = 0; i < NP; ++i)
            m_particles[i](di) += sigma * normal(m_rng);
    }
};

// Factory function since NP cannot be deduced via CTAD
template <std::size_t NP, typename Dynamics, typename Measurement, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Rng = std::mt19937_64>
auto make_particle_filter(Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{})
{
    return particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement, systematic_resampling, Rng>(std::move(dynamics), std::move(measurement), std::move(config), std::move(rng));
}

// Static assert helpers
namespace detail
{

struct pf_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>&, const Vector<double, 1>&) const { return Vector<double, 2>::Zero(); }
};

struct pf_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>&) const { return Vector<double, 1>::Zero(); }
};

using pf_test_type = particle_filter<double, 2, 1, 1, 10, pf_sa_dynamics, pf_sa_measurement>;

}

static_assert(ObserverPolicy<detail::pf_test_type>);
static_assert(!CovarianceObserver<detail::pf_test_type>);

}

#endif
