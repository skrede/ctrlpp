#ifndef HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H
#define HPP_GUARD_CTRLPP_DSP_DISCRETE_FILTER_H

/// @brief Discrete filter concept: minimal interface for composable digital filters.
///
/// @cite oppenheim1997 -- Oppenheim & Willsky, "Signals and Systems", 1997

#include <concepts>

namespace ctrlpp {

template<typename F>
concept discrete_filter = requires(F f, typename F::scalar_type x) {
    { f.process(x) } -> std::convertible_to<typename F::scalar_type>;
};

}

#endif
