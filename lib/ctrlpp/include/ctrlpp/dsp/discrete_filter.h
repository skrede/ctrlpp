#ifndef HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H
#define HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H

// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Secs. 2.3 and 2.5, Eqs. (2.49) and (2.84) (convolution sum, difference equations)

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
