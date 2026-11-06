#ifndef HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_H

// @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/particle_filter_types.h"

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"

#include "ctrlpp/estimation/resampling/resampling_strategy.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"

#include "ctrlpp/estimation/detail/pf_weights.h"

#include <array>
#include <random>
#include <cstddef>
#include <utility>
#include <concepts>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement,
        typename Resampler = systematic_resampling, typename Rng = std::mt19937_64, typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY> && std::uniform_random_bit_generator<Rng> &&
        resampling_strategy<Resampler, Rng, NP> && pf_likelihood_model<Likelihood, Scalar, NY>
class particle_filter : private detail::pf_weights<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    static_assert(NP > 0, "Particle count NP must be positive");

    using base = detail::pf_weights<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>;
    using typename base::sampling_factors;

    // Only the filter can make the passkey, so create is the only way to originate a filter.
    class validated_tag
    {
        friend particle_filter;
        validated_tag() = default;
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
        : base{std::move(dynamics), std::move(measurement), std::move(config), std::move(factors), std::move(resampler), std::move(rng), std::move(likelihood)}
    {
    }

    static ctrlpp::expected<particle_filter, filter_error> create(
            Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
    {
        return create(std::move(dynamics), std::move(measurement), std::move(config), Resampler{}, std::move(rng), std::move(likelihood));
    }

    static ctrlpp::expected<particle_filter, filter_error> create(
            Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config, Resampler resampler, Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
    {
        auto factors = base::validated_sampling_factors(config);
        if(!factors)
            return ctrlpp::unexpected(factors.error());
        return ctrlpp::expected<particle_filter, filter_error>{std::in_place, validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config), std::move(*factors),
                std::move(resampler), std::move(rng), std::move(likelihood)};
    }

    // @cite gordon1993 -- Gordon et al., "Novel approach to nonlinear/non-Gaussian Bayesian state estimation", 1993, Sec. 3.1
    void predict(const input_vector_t& u) { this->propagate_particles(u); }

    // @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002, Sec. V-A, Alg. 3, Eq. 63
    ctrlpp::expected<void, pf_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));
        if(this->m_weight_mode == weight_representation::log)
            return this->update_log(z);
        return this->update_linear(z);
    }

    state_vector_t state() const
    {
        if(this->m_extraction == extraction_method::map)
            return this->map_estimate();
        return this->weighted_mean();
    }

    using base::weighted_mean;
    using base::covariance;
    using base::map_estimate;

    const std::array<state_vector_t, NP>& particles() const { return this->m_particles; }

    pf_health health() const { return this->m_health; }

private:
    // Classifies without touching a member: the carried particles first, the measurement last.
    ctrlpp::expected<void, pf_update_error> check_step(const output_vector_t& z) const
    {
        for(const auto& particle : this->m_particles)
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
            this->m_health = pf_health::non_finite_estimate;
        return fault;
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
