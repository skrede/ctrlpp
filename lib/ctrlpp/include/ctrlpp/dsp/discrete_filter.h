#ifndef HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H
#define HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H

#include <concepts>

namespace ctrlpp {

template<typename F>
concept discrete_filter = requires(F f, typename F::scalar_type x) {
    { f.process(x) } -> std::convertible_to<typename F::scalar_type>;
};

}

#endif
