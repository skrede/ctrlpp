#ifndef HPP_GUARD_CTRLPP_UTIL_CONCEPTS_H
#define HPP_GUARD_CTRLPP_UTIL_CONCEPTS_H

#include <concepts>

namespace ctrlpp
{

template <typename S>
concept ctrlpp_floating_scalar = std::floating_point<S>;

}

#endif
