#ifndef HPP_GUARD_CTRLPP_ESTIMATION_RESAMPLING_SYSTEMATIC_RESAMPLING_H
#define HPP_GUARD_CTRLPP_ESTIMATION_RESAMPLING_SYSTEMATIC_RESAMPLING_H

/// @brief Systematic resampling for particle filters.
///
/// @cite arulampalam2002 -- Arulampalam et al., "A Tutorial on Particle Filters", 2002

#include <array>
#include <random>
#include <cstddef>

namespace ctrlpp {

struct systematic_resampling
{
    template <typename Scalar, std::size_t NP, typename Rng>
    void resample(const std::array<Scalar, NP> &weights, std::array<std::size_t, NP> &indices, Rng &rng) const
    {
        std::array<Scalar, NP> cumsum;
        cumsum[0] = weights[0];
        for(std::size_t i = 1; i < NP; ++i)
            cumsum[i] = cumsum[i - 1] + weights[i];

        std::uniform_real_distribution<Scalar> dist(Scalar{0}, Scalar{1} / static_cast<Scalar>(NP));
        Scalar u = dist(rng);

        std::size_t j = 0;
        for(std::size_t i = 0; i < NP; ++i)
        {
            Scalar threshold = u + static_cast<Scalar>(i) / static_cast<Scalar>(NP);
            while(j < NP - 1 && cumsum[j] < threshold)
                ++j;
            indices[i] = j;
        }
    }
};

}

#endif
