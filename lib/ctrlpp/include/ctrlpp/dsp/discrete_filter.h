#ifndef HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H
#define HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H

/// @brief Discrete filter concept: minimal interface for composable digital filters.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 2 (LTI difference equations and convolution)

#include <concepts>

namespace ctrlpp
{

template <typename F>
concept discrete_filter = requires(F f, typename F::scalar_type x) {
    { f.process(x) } -> std::convertible_to<typename F::scalar_type>;
};

template <typename F, typename Vec>
concept vector_discrete_filter = requires(F f, Vec v) {
    { f.process(v) } -> std::convertible_to<Vec>;
};

}

#endif
