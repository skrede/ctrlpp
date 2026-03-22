#ifndef HPP_GUARD_CTRLPP_RESAMPLING_MULTINOMIAL_RESAMPLING_H
#define HPP_GUARD_CTRLPP_RESAMPLING_MULTINOMIAL_RESAMPLING_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>

namespace ctrlpp {

struct multinomial_resampling {
    template<typename Scalar, std::size_t NP, typename Rng>
    void resample(const std::array<Scalar, NP>& weights,
                  std::array<std::size_t, NP>& indices,
                  Rng& rng) const
    {
        std::array<Scalar, NP> cumsum;
        cumsum[0] = weights[0];
        for (std::size_t i = 1; i < NP; ++i) {
            cumsum[i] = cumsum[i - 1] + weights[i];
        }

        std::array<Scalar, NP> samples;
        std::uniform_real_distribution<Scalar> dist(Scalar{0}, Scalar{1});
        for (std::size_t i = 0; i < NP; ++i) {
            samples[i] = dist(rng);
        }
        std::sort(samples.begin(), samples.end());

        std::size_t j = 0;
        for (std::size_t i = 0; i < NP; ++i) {
            while (j < NP - 1 && cumsum[j] < samples[i]) {
                ++j;
            }
            indices[i] = j;
        }
    }
};

}

#endif
