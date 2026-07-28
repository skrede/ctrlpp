#ifndef HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H

/// @brief Bootstrap SIR particle filter with ESS-adaptive resampling and roughening.
///
/// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/observer_policy.h"
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
#include <concepts>
#include <algorithm>

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

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct pf_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar ess_threshold{Scalar{-1}};
    Scalar roughening_scale{Scalar{0.2}};
    extraction_method extraction{extraction_method::weighted_mean};
    weight_representation weights{weight_representation::log};
};

namespace detail
{

template <ctrlpp_floating_scalar Scalar, std::size_t NY>
struct gaussian_likelihood
{
    /// Log Gaussian likelihood of the innovation z - z_pred under measurement
    /// covariance R, given precomputed R^-1 and log((2 pi)^NY det R).
    Scalar operator()(const Vector<Scalar, NY>& z, const Vector<Scalar, NY>& z_pred,
                      const Matrix<Scalar, NY, NY>& R_inv, Scalar log_det_2piR) const
    {
        Vector<Scalar, NY> innov = z - z_pred;
        Scalar mahal = (innov.transpose() * R_inv * innov)(0, 0);
        return Scalar{-0.5} * mahal - Scalar{0.5} * log_det_2piR;
    }
};

}

/// @brief Constrains a measurement-likelihood policy usable by particle_filter.
///
/// A policy maps an innovation (z, z_pred) plus the precomputed measurement
/// precision R^-1 and log normalizer to a scalar log-likelihood.
template <typename L, typename Scalar, std::size_t NY>
concept pf_likelihood_model = requires(const L& l, const Vector<Scalar, NY>& v, const Matrix<Scalar, NY, NY>& Rinv, Scalar s) {
    { l(v, v, Rinv, s) } -> std::convertible_to<Scalar>;
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement, typename Resampler = systematic_resampling, typename Rng = std::mt19937_64, typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && std::uniform_random_bit_generator<Rng> && resampling_strategy<Resampler, Rng, NP> && pf_likelihood_model<Likelihood, Scalar, NY>
class particle_filter
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    static_assert(NP > 0, "Particle count NP must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct pf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;

    particle_filter(Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_resampler{}
        , m_rng{std::move(rng)}
        , m_likelihood{std::move(likelihood)}
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

    particle_filter(Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Resampler resampler, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_resampler{std::move(resampler)}
        , m_rng{std::move(rng)}
        , m_likelihood{std::move(likelihood)}
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

    auto weighted_mean() const -> state_vector_t
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

    /// @brief Weighted second central moment of the particle set: the posterior
    /// covariance the cloud and its weights represent.
    ///
    /// This is the uncertainty half of the estimate, and without it a caller has
    /// no way to learn how much to trust `state()`. The particle array alone
    /// does not answer the question: its unweighted dispersion ignores the
    /// weights entirely and therefore describes the PRIOR spread whenever the
    /// last update did not trigger a resampling.
    ///
    /// Three properties of the definition, each chosen rather than defaulted:
    ///
    ///  * The centre is the WEIGHTED MEAN, always, including when the extraction
    ///    method is the maximum a posteriori particle. A second moment about any
    ///    other point is larger than the covariance and is not one; a caller who
    ///    wants the dispersion about the MAP estimate can form it from
    ///    `particles()` and `map_estimate()`.
    ///  * There is no Bessel correction. The weights sum to one, so this is the
    ///    weighted second moment, matching the convention the unscented filter's
    ///    own sigma-point covariance uses.
    ///  * Uniform weights need no special case and get none. The expression then
    ///    reduces exactly to the plain second moment of the particles about their
    ///    plain mean, which is the honest answer after the weight-degeneracy
    ///    guard fires: the measurement carried no information, so the reported
    ///    uncertainty is the dispersion the filter was already carrying.
    ///
    /// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. III-A
    auto covariance() const -> Matrix<Scalar, NX, NX>
    {
        const state_vector_t mean = weighted_mean();
        Matrix<Scalar, NX, NX> P = Matrix<Scalar, NX, NX>::Zero();

        if(m_weight_mode == weight_representation::log)
        {
            auto lin = log_to_linear();
            for(std::size_t i = 0; i < NP; ++i)
            {
                auto diff = (m_particles[i] - mean).eval();
                P += lin[i] * diff * diff.transpose();
            }
        }
        else
        {
            for(std::size_t i = 0; i < NP; ++i)
            {
                auto diff = (m_particles[i] - mean).eval();
                P += m_linear_weights[i] * diff * diff.transpose();
            }
        }

        return detail::symmetrize(P);
    }

    auto map_estimate() const -> state_vector_t
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

    auto particles() const -> const std::array<state_vector_t, NP>& { return m_particles; }

private:
    Dynamics m_dynamics;
    Measurement m_measurement;
    Resampler m_resampler;
    Rng m_rng;
    Likelihood m_likelihood;

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

    /// @brief Compute measurement log-likelihood of the innovation via the policy.
    ///
    /// The default policy reproduces the Gaussian log-likelihood exactly; a
    /// custom policy (for example a bearing wrap) sees the same precomputed
    /// precision m_R_inv and log normalizer m_log_det_2piR.
    ///
    /// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 63
    Scalar log_likelihood(const output_vector_t& z, const output_vector_t& z_pred) const
    {
        return m_likelihood(z, z_pred, m_R_inv, m_log_det_2piR);
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
    /// When all particles have negligible likelihood (max log-weight is -inf),
    /// resets to uniform weights to recover from particle depletion.
    ///
    /// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. III-A
    void normalize_log_weights()
    {
        Scalar max_log_w = *std::max_element(m_log_weights.begin(), m_log_weights.end());
        if(!std::isfinite(max_log_w))
        {
            Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
            m_log_weights.fill(log_uniform);
            return;
        }
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
    ///
    /// On total underflow (every likelihood collapsed to zero, so the sum is
    /// not positive) the weights are reset to uniform 1/NP, mirroring the log
    /// path's reset, so the filter recovers to the plain particle mean rather
    /// than leaving stale or zero weights.
    void normalize_linear_weights()
    {
        Scalar sum = Scalar{0};
        for(auto w : m_linear_weights)
            sum += w;
        if(sum > Scalar{0})
            for(auto& w : m_linear_weights)
                w /= sum;
        else
            m_linear_weights.fill(Scalar{1} / static_cast<Scalar>(NP));
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
template <std::size_t NP, typename Dynamics, typename Measurement, ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Rng = std::mt19937_64, typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
auto make_particle_filter(Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
{
    return particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement, systematic_resampling, Rng, Likelihood>(std::move(dynamics), std::move(measurement), std::move(config), std::move(rng), std::move(likelihood));
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
