#ifndef HPP_GUARD_CTRLPP_RESAMPLING_RESAMPLING_STRATEGY_H
#define HPP_GUARD_CTRLPP_RESAMPLING_RESAMPLING_STRATEGY_H

#include <array>
#include <cstddef>

namespace ctrlpp {

template <typename R, typename Rng, std::size_t NP>
concept resampling_strategy = requires(const R &r, const std::array<double, NP> &weights, std::array<std::size_t, NP> &indices, Rng &rng)
{
    r.resample(weights, indices, rng);
};

}

#endif
