#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_PF_PARTICLES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_PF_PARTICLES_H

#include "ctrlpp/estimation/detail/pf_storage.h"

#include "ctrlpp/estimation/particle_filter_types.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <array>
#include <cmath>
#include <random>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement, typename Resampler, typename Rng,
        typename Likelihood>
class pf_particles : protected pf_storage<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>
{
    using base = pf_storage<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>;

public:
    using typename base::state_vector_t;
    using typename base::input_vector_t;
    using typename base::output_vector_t;

    state_vector_t weighted_mean() const
    {
        state_vector_t mean = state_vector_t::Zero();
        if(this->m_weight_mode == weight_representation::log)
        {
            auto lin = log_to_linear();
            for(std::size_t i = 0; i < NP; ++i)
                mean += lin[i] * this->m_particles[i];
        }
        else
        {
            for(std::size_t i = 0; i < NP; ++i)
                mean += this->m_linear_weights[i] * this->m_particles[i];
        }
        return mean;
    }

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. III-A
    Matrix<Scalar, NX, NX> covariance() const
    {
        const state_vector_t mean = weighted_mean();
        Matrix<Scalar, NX, NX> P = Matrix<Scalar, NX, NX>::Zero();

        if(this->m_weight_mode == weight_representation::log)
        {
            auto lin = log_to_linear();
            for(std::size_t i = 0; i < NP; ++i)
            {
                auto diff = (this->m_particles[i] - mean).eval();
                P += lin[i] * diff * diff.transpose();
            }
        }
        else
        {
            for(std::size_t i = 0; i < NP; ++i)
            {
                auto diff = (this->m_particles[i] - mean).eval();
                P += this->m_linear_weights[i] * diff * diff.transpose();
            }
        }

        return detail::symmetrize(P);
    }

    state_vector_t map_estimate() const
    {
        std::size_t best = 0;
        if(this->m_weight_mode == weight_representation::log)
        {
            for(std::size_t i = 1; i < NP; ++i)
            {
                if(this->m_log_weights[i] > this->m_log_weights[best])
                    best = i;
            }
        }
        else
        {
            for(std::size_t i = 1; i < NP; ++i)
            {
                if(this->m_linear_weights[i] > this->m_linear_weights[best])
                    best = i;
            }
        }
        return this->m_particles[best];
    }

protected:
    using base::base;

    std::array<Scalar, NP> log_to_linear() const
    {
        std::array<Scalar, NP> lin;
        for(std::size_t i = 0; i < NP; ++i)
            lin[i] = std::exp(this->m_log_weights[i]);
        return lin;
    }

    void note_resampling_rolled_back()
    {
        if(this->m_health < pf_health::resampling_rolled_back)
            this->m_health = pf_health::resampling_rolled_back;
    }

    void reindex_particles(const std::array<std::size_t, NP>& indices)
    {
        this->m_previous_particles = this->m_particles;
        for(std::size_t i = 0; i < NP; ++i)
            this->m_particles[i] = this->m_previous_particles[indices[i]];
    }

    // On a non-finite roughening the set before the resampling is restored, and the caller keeps
    // the step's weights, which are the posterior over that set.
    bool roughened_particles_kept()
    {
        for(const auto& particle : this->m_particles)
        {
            if(!particle.allFinite())
            {
                this->m_particles = this->m_previous_particles;
                note_resampling_rolled_back();
                return false;
            }
        }
        return true;
    }

    // @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Sec. 5
    void apply_roughening()
    {
        if(this->m_roughening_scale <= Scalar{0})
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
        Scalar min_val = this->m_particles[0](di);
        Scalar max_val = min_val;
        for(std::size_t i = 1; i < NP; ++i)
        {
            Scalar v = this->m_particles[i](di);
            min_val = std::min(min_val, v);
            max_val = std::max(max_val, v);
        }
        Scalar sigma = this->m_roughening_scale * (max_val - min_val) * np_factor;

        // A collapsed coordinate falls back to the process-noise scale.
        if(sigma <= Scalar{0})
            sigma = this->m_roughening_scale * std::sqrt(this->m_Q(di, di));
        return sigma;
    }

    void apply_particle_jitter(std::size_t d, Scalar sigma, std::normal_distribution<Scalar>& normal)
    {
        const auto di = static_cast<Eigen::Index>(d);
        for(std::size_t i = 0; i < NP; ++i)
            this->m_particles[i](di) += sigma * normal(this->m_rng);
    }
};

}

#endif
