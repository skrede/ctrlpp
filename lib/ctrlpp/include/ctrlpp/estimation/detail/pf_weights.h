#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_PF_WEIGHTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_PF_WEIGHTS_H

#include "ctrlpp/estimation/detail/pf_particles.h"

#include "ctrlpp/estimation/particle_filter_types.h"

#include "ctrlpp/expected.h"

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement, typename Resampler, typename Rng,
        typename Likelihood>
class pf_weights : protected pf_particles<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>
{
    using base = pf_particles<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>;

public:
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;

protected:
    using base::base;

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 63
    Scalar log_likelihood(const output_vector_t& z, const output_vector_t& z_pred) const { return this->m_likelihood(z, z_pred, this->m_R_inv, this->m_log_det_2piR); }

    void note_weights_reset()
    {
        if(this->m_health == pf_health::ok)
            this->m_health = pf_health::weights_reset;
    }

    void resample_log()
    {
        auto lin = this->log_to_linear();
        std::array<std::size_t, NP> indices{};
        this->m_resampler.resample(lin, indices, this->m_rng);
        this->reindex_particles(indices);
        this->apply_roughening();
        if(!this->roughened_particles_kept())
            return;

        Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
        this->m_log_weights.fill(log_uniform);
    }

    void resample_linear()
    {
        std::array<std::size_t, NP> indices{};
        this->m_resampler.resample(this->m_linear_weights, indices, this->m_rng);
        this->reindex_particles(indices);
        this->apply_roughening();
        if(!this->roughened_particles_kept())
            return;

        Scalar uniform = Scalar{1} / static_cast<Scalar>(NP);
        this->m_linear_weights.fill(uniform);
    }

    // A -inf log-weight is a legitimate zero likelihood, so only NaN and +inf refuse the step.
    bool stage_log_weights(const output_vector_t& z)
    {
        bool admissible = true;
        for(std::size_t i = 0; i < NP; ++i)
        {
            output_vector_t z_pred = this->m_measurement(this->m_particles[i]);
            this->m_staged_weights[i] = this->m_log_weights[i] + log_likelihood(z, z_pred);
            admissible = admissible && !std::isnan(this->m_staged_weights[i]) && this->m_staged_weights[i] != std::numeric_limits<Scalar>::infinity();
        }
        return admissible;
    }

    // Log-sum-exp normalization.
    void normalize_log_weights()
    {
        Scalar max_log_w = *std::max_element(this->m_log_weights.begin(), this->m_log_weights.end());
        if(!std::isfinite(max_log_w))
        {
            Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
            this->m_log_weights.fill(log_uniform);
            note_weights_reset();
            return;
        }
        Scalar sum_exp = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
            sum_exp += std::exp(this->m_log_weights[i] - max_log_w);
        Scalar log_sum = max_log_w + std::log(sum_exp);
        for(auto& lw : this->m_log_weights)
            lw -= log_sum;
    }

    // The total is summed in the order normalization sums it, so the two agree on finiteness.
    bool stage_linear_weights(const output_vector_t& z)
    {
        Scalar total = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
        {
            output_vector_t z_pred = this->m_measurement(this->m_particles[i]);
            Scalar ll = log_likelihood(z, z_pred);
            this->m_staged_weights[i] = this->m_linear_weights[i] * std::exp(ll);
            total += this->m_staged_weights[i];
        }
        return std::isfinite(total);
    }

    void normalize_linear_weights()
    {
        Scalar sum = Scalar{0};
        for(auto w : this->m_linear_weights)
            sum += w;
        if(sum > Scalar{0})
        {
            for(auto& w : this->m_linear_weights)
                w /= sum;
            return;
        }
        this->m_linear_weights.fill(Scalar{1} / static_cast<Scalar>(NP));
        note_weights_reset();
    }

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Eq. 51
    Scalar compute_ess_from_log() const
    {
        Scalar sum_w2 = Scalar{0};
        for(std::size_t i = 0; i < NP; ++i)
        {
            Scalar w = std::exp(this->m_log_weights[i]);
            sum_w2 += w * w;
        }
        return Scalar{1} / sum_w2;
    }

    Scalar compute_ess_from_linear() const
    {
        Scalar sum_w2 = Scalar{0};
        for(auto w : this->m_linear_weights)
            sum_w2 += w * w;
        return Scalar{1} / sum_w2;
    }

    ctrlpp::expected<void, pf_update_error> update_log(const output_vector_t& z)
    {
        if(!stage_log_weights(z))
            return ctrlpp::unexpected(pf_update_error::non_finite_result);
        this->m_log_weights = this->m_staged_weights;
        normalize_log_weights();
        if(compute_ess_from_log() < this->m_ess_threshold)
            resample_log();
        return {};
    }

    ctrlpp::expected<void, pf_update_error> update_linear(const output_vector_t& z)
    {
        if(!stage_linear_weights(z))
            return ctrlpp::unexpected(pf_update_error::non_finite_result);
        this->m_linear_weights = this->m_staged_weights;
        normalize_linear_weights();
        if(compute_ess_from_linear() < this->m_ess_threshold)
            resample_linear();
        return {};
    }
};

}

#endif
