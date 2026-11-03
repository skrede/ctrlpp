#ifndef HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H

// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/resampling/resampling_strategy.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"

#include "ctrlpp/estimation/detail/filter_configuration.h"

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

enum class pf_update_error
{
    non_finite_state,
    non_finite_measurement,
    non_finite_result
};

// Ordered by severity; the status latches, so it never returns to a lower rung.
enum class pf_health
{
    ok,
    weights_reset,
    resampling_rolled_back,
    non_finite_estimate
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
    Scalar roughening_scale{static_cast<Scalar>(0.2)};
    extraction_method extraction{extraction_method::weighted_mean};
    weight_representation weights{weight_representation::log};
};

namespace detail
{

template <ctrlpp_floating_scalar Scalar, std::size_t NY>
struct gaussian_likelihood
{
    Scalar operator()(const Vector<Scalar, NY>& z, const Vector<Scalar, NY>& z_pred, const Matrix<Scalar, NY, NY>& R_inv, Scalar log_det_2piR) const
    {
        Vector<Scalar, NY> innov = z - z_pred;
        Scalar mahal = (innov.transpose() * R_inv * innov)(0, 0);
        return Scalar{-0.5} * mahal - Scalar{0.5} * log_det_2piR;
    }
};

}

template <typename L, typename Scalar, std::size_t NY>
concept pf_likelihood_model = requires(const L& l, const Vector<Scalar, NY>& v, const Matrix<Scalar, NY, NY>& Rinv, Scalar s) {
    { l(v, v, Rinv, s) } -> std::convertible_to<Scalar>;
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement,
        typename Resampler = systematic_resampling, typename Rng = std::mt19937_64, typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && std::uniform_random_bit_generator<Rng> &&
        resampling_strategy<Resampler, Rng, NP> && pf_likelihood_model<Likelihood, Scalar, NY>
class particle_filter
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    static_assert(NP > 0, "Particle count NP must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);

    // Only the filter can make the passkey, so create is the only way to originate a filter.
    class validated_tag
    {
        friend particle_filter;
        validated_tag() = default;
    };

    // Factors S with S S^T equal to Q (process) and to P0 (prior).
    struct sampling_factors
    {
        Matrix<Scalar, NX, NX> process;
        Matrix<Scalar, NX, NX> prior;
    };

public:
    using observer_tag = struct pf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;

    // Public so ctrlpp::expected can construct the filter in the caller's
    // storage; unreachable without the validated_tag only the filter can make.
    particle_filter(validated_tag, Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, sampling_factors factors, Resampler resampler, Rng rng,
            Likelihood likelihood)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_resampler{std::move(resampler)}
        , m_rng{std::move(rng)}
        , m_likelihood{std::move(likelihood)}
        , m_Q{std::move(config.Q)}
        , m_Q_sqrt{std::move(factors.process)}
        , m_R{std::move(config.R)}
        , m_R_inv{m_R.colPivHouseholderQr().inverse()}
        , m_log_det_2piR{compute_log_det_2piR()}
        , m_ess_threshold{config.ess_threshold < Scalar{0} ? static_cast<Scalar>(NP) / Scalar{2} : config.ess_threshold}
        , m_roughening_scale{config.roughening_scale}
        , m_extraction{config.extraction}
        , m_weight_mode{config.weights}
        , m_health{pf_health::ok}
        , m_previous_particles{}
        , m_log_weights{}
        , m_linear_weights{}
        , m_staged_weights{}
    {
        initialize_particles(config.x0, factors.prior);
    }

    static ctrlpp::expected<particle_filter, filter_error> create(
            Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
    {
        return create(std::move(dynamics), std::move(measurement), std::move(config), Resampler{}, std::move(rng), std::move(likelihood));
    }

    static ctrlpp::expected<particle_filter, filter_error> create(
            Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Resampler resampler, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
    {
        auto factors = validated_sampling_factors(config);
        if(!factors)
            return ctrlpp::unexpected(factors.error());
        return ctrlpp::expected<particle_filter, filter_error>{std::in_place, validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config), std::move(*factors),
                std::move(resampler), std::move(rng), std::move(likelihood)};
    }

    // @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Eq. 6
    void predict(const input_vector_t& u) { propagate_particles(u); }

    // @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Alg. 1
    ctrlpp::expected<void, pf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));
        if(m_weight_mode == weight_representation::log)
            return update_log(z);
        return update_linear(z);
    }

    state_vector_t state() const
    {
        if(m_extraction == extraction_method::map)
            return map_estimate();
        return weighted_mean();
    }

    state_vector_t weighted_mean() const
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

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. III-A
    Matrix<Scalar, NX, NX> covariance() const
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

    state_vector_t map_estimate() const
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

    const std::array<state_vector_t, NP>& particles() const { return m_particles; }

    pf_health health() const { return m_health; }

private:
    Dynamics m_dynamics;
    Measurement m_measurement;
    Resampler m_resampler;
    Rng m_rng;
    Likelihood m_likelihood;

    Matrix<Scalar, NX, NX> m_Q;
    Matrix<Scalar, NX, NX> m_Q_sqrt;
    Matrix<Scalar, NY, NY> m_R;
    Matrix<Scalar, NY, NY> m_R_inv;
    Scalar m_log_det_2piR;
    Scalar m_ess_threshold;
    Scalar m_roughening_scale;
    extraction_method m_extraction;
    weight_representation m_weight_mode;
    pf_health m_health;

    std::array<state_vector_t, NP> m_particles;
    // The set before the current resampling: the resampling gathers from it,
    // and a roughening that leaves a coordinate non-finite is undone from it.
    std::array<state_vector_t, NP> m_previous_particles;
    std::array<Scalar, NP> m_log_weights;
    std::array<Scalar, NP> m_linear_weights;
    // Staged so a refused update leaves the weights untouched; a member rather
    // than a local keeps the update's stack frame independent of NP.
    std::array<Scalar, NP> m_staged_weights;

    static ctrlpp::expected<sampling_factors, filter_error> validated_sampling_factors(pf_config<Scalar, NX, NU, NY>& config)
    {
        const auto Q = detail::checked_covariance(config.Q, detail::process_noise_faults);
        if(!Q)
            return ctrlpp::unexpected(Q.error());
        const auto R = detail::checked_covariance(config.R, detail::measurement_noise_faults);
        if(!R)
            return ctrlpp::unexpected(R.error());
        if(!detail::nonsingular_spectrum(*R))
            return ctrlpp::unexpected(filter_error::singular_measurement_noise);
        if(!config.x0.allFinite())
            return ctrlpp::unexpected(filter_error::non_finite_initial_state);
        const auto P0 = detail::checked_covariance(config.P0, detail::initial_covariance_faults);
        if(!P0)
            return ctrlpp::unexpected(P0.error());
        return sampling_factors{sampling_factor(config.Q, *Q), sampling_factor(config.P0, *P0)};
    }

    // Eigen's LLT leaves the input's entries in the columns after a non-positive pivot, so its
    // factor is read only after a successful factorization.
    static Matrix<Scalar, NX, NX> sampling_factor(const Matrix<Scalar, NX, NX>& M, const detail::psd_decomposition_result<Scalar, nx>& decomposition)
    {
        Eigen::LLT<Matrix<Scalar, NX, NX>> llt(M);
        if(llt.info() == Eigen::Success)
            return llt.matrixL();
        return detail::psd_factor(decomposition);
    }

    Scalar compute_log_det_2piR() const
    {
        Scalar log_det = m_R.colPivHouseholderQr().logAbsDeterminant();
        return static_cast<Scalar>(NY) * std::log(Scalar{2} * std::numbers::pi_v<Scalar>) + log_det;
    }

    void initialize_particles(const state_vector_t& x0, const Matrix<Scalar, NX, NX>& P0_sqrt)
    {
        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});

        for(std::size_t i = 0; i < NP; ++i)
        {
            state_vector_t noise;
            for(std::size_t d = 0; d < NX; ++d)
                noise(static_cast<Eigen::Index>(d)) = normal(m_rng);
            m_particles[i] = x0 + P0_sqrt * noise;
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
            noise(static_cast<Eigen::Index>(d)) = normal(m_rng);
        return m_Q_sqrt * noise;
    }

    void propagate_particles(const input_vector_t& u)
    {
        for(std::size_t i = 0; i < NP; ++i)
        {
            m_particles[i] = m_dynamics(m_particles[i], u);
            m_particles[i] += sample_process_noise();
        }
    }

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 63
    Scalar log_likelihood(const output_vector_t& z, const output_vector_t& z_pred) const { return m_likelihood(z, z_pred, m_R_inv, m_log_det_2piR); }

    // Classifies without touching a member: the carried particles first, the measurement last.
    ctrlpp::expected<void, pf_update_error> check_step(const output_vector_t& z) const
    {
        for(const auto& particle : m_particles)
        {
            if(!particle.allFinite())
                return ctrlpp::unexpected(pf_update_error::non_finite_state);
        }
        if(!z.allFinite())
            return ctrlpp::unexpected(pf_update_error::non_finite_measurement);
        return {};
    }

    pf_update_error latch_health(pf_update_error fault)
    {
        if(fault != pf_update_error::non_finite_measurement)
            m_health = pf_health::non_finite_estimate;
        return fault;
    }

    void note_weights_reset()
    {
        if(m_health == pf_health::ok)
            m_health = pf_health::weights_reset;
    }

    void note_resampling_rolled_back()
    {
        if(m_health < pf_health::resampling_rolled_back)
            m_health = pf_health::resampling_rolled_back;
    }

    // A -inf log-weight is a legitimate zero likelihood, so only NaN and +inf refuse the step.
    bool stage_log_weights(const output_vector_t& z)
    {
        bool admissible = true;
        for(std::size_t i = 0; i < NP; ++i)
        {
            output_vector_t z_pred = m_measurement(m_particles[i]);
            m_staged_weights[i] = m_log_weights[i] + log_likelihood(z, z_pred);
            admissible = admissible && !std::isnan(m_staged_weights[i]) && m_staged_weights[i] != std::numeric_limits<Scalar>::infinity();
        }
        return admissible;
    }

    // Log-sum-exp normalization.
    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. III-A
    void normalize_log_weights()
    {
        Scalar max_log_w = *std::max_element(m_log_weights.begin(), m_log_weights.end());
        if(!std::isfinite(max_log_w))
        {
            Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
            m_log_weights.fill(log_uniform);
            note_weights_reset();
            return;
        }
        Scalar sum_exp = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
            sum_exp += std::exp(m_log_weights[i] - max_log_w);
        Scalar log_sum = max_log_w + std::log(sum_exp);
        for(auto& lw : m_log_weights)
            lw -= log_sum;
    }

    // The total is summed in the order normalization sums it, so the two agree on finiteness.
    bool stage_linear_weights(const output_vector_t& z)
    {
        Scalar total = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
        {
            output_vector_t z_pred = m_measurement(m_particles[i]);
            Scalar ll = log_likelihood(z, z_pred);
            m_staged_weights[i] = m_linear_weights[i] * std::exp(ll);
            total += m_staged_weights[i];
        }
        return std::isfinite(total);
    }

    void normalize_linear_weights()
    {
        Scalar sum = Scalar{0};
        for(auto w : m_linear_weights)
            sum += w;
        if(sum > Scalar{0})
        {
            for(auto& w : m_linear_weights)
                w /= sum;
            return;
        }
        m_linear_weights.fill(Scalar{1} / static_cast<Scalar>(NP));
        note_weights_reset();
    }

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 51
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

    Scalar compute_ess_from_linear() const
    {
        Scalar sum_w2 = Scalar{0};
        for(auto w : m_linear_weights)
            sum_w2 += w * w;
        return Scalar{1} / sum_w2;
    }

    void resample_log()
    {
        auto lin = log_to_linear();
        std::array<std::size_t, NP> indices{};
        m_resampler.resample(lin, indices, m_rng);
        reindex_particles(indices);
        apply_roughening();
        if(!roughened_particles_kept())
            return;

        Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
        m_log_weights.fill(log_uniform);
    }

    void resample_linear()
    {
        std::array<std::size_t, NP> indices{};
        m_resampler.resample(m_linear_weights, indices, m_rng);
        reindex_particles(indices);
        apply_roughening();
        if(!roughened_particles_kept())
            return;

        Scalar uniform = Scalar{1} / static_cast<Scalar>(NP);
        m_linear_weights.fill(uniform);
    }

    ctrlpp::expected<void, pf_update_error> update_log(const output_vector_t& z)
    {
        if(!stage_log_weights(z))
            return ctrlpp::unexpected(pf_update_error::non_finite_result);
        m_log_weights = m_staged_weights;
        normalize_log_weights();
        if(compute_ess_from_log() < m_ess_threshold)
            resample_log();
        return {};
    }

    ctrlpp::expected<void, pf_update_error> update_linear(const output_vector_t& z)
    {
        if(!stage_linear_weights(z))
            return ctrlpp::unexpected(pf_update_error::non_finite_result);
        m_linear_weights = m_staged_weights;
        normalize_linear_weights();
        if(compute_ess_from_linear() < m_ess_threshold)
            resample_linear();
        return {};
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
        m_previous_particles = m_particles;
        for(std::size_t i = 0; i < NP; ++i)
            m_particles[i] = m_previous_particles[indices[i]];
    }

    // On a non-finite roughening the set before the resampling is restored, and the caller keeps
    // the step's weights, which are the posterior over that set.
    bool roughened_particles_kept()
    {
        for(const auto& particle : m_particles)
        {
            if(!particle.allFinite())
            {
                m_particles = m_previous_particles;
                note_resampling_rolled_back();
                return false;
            }
        }
        return true;
    }

    // @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Sec. 5
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
        const auto di = static_cast<Eigen::Index>(d);
        Scalar min_val = m_particles[0](di);
        Scalar max_val = min_val;
        for(std::size_t i = 1; i < NP; ++i)
        {
            Scalar v = m_particles[i](di);
            min_val = std::min(min_val, v);
            max_val = std::max(max_val, v);
        }
        Scalar sigma = m_roughening_scale * (max_val - min_val) * np_factor;

        // A collapsed coordinate falls back to the process-noise scale.
        if(sigma <= Scalar{0})
            sigma = m_roughening_scale * std::sqrt(m_Q(di, di));
        return sigma;
    }

    void apply_particle_jitter(std::size_t d, Scalar sigma, std::normal_distribution<Scalar>& normal)
    {
        const auto di = static_cast<Eigen::Index>(d);
        for(std::size_t i = 0; i < NP; ++i)
            m_particles[i](di) += sigma * normal(m_rng);
    }
};

// NP leads because it cannot be deduced from the arguments.
template <std::size_t NP, typename Dynamics, typename Measurement, ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Rng = std::mt19937_64,
        typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
ctrlpp::expected<particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement, systematic_resampling, Rng, Likelihood>, filter_error> make_particle_filter(
        Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
{
    using filter_type = particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement, systematic_resampling, Rng, Likelihood>;
    return filter_type::create(std::move(dynamics), std::move(measurement), std::move(config), std::move(rng), std::move(likelihood));
}

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
