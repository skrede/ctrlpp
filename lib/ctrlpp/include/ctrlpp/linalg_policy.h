#ifndef HPP_GUARD_CTRLPP_LINALG_POLICY_H
#define HPP_GUARD_CTRLPP_LINALG_POLICY_H

#include <cstddef>
#include <concepts>
#include <type_traits>

namespace ctrlpp {

template<typename P>
concept LinalgPolicy = requires {
    typename P::template matrix_type<double, 2, 2>;
    typename P::template matrix_type<double, 2, 3>;
    typename P::template matrix_type<double, 3, 2>;
    typename P::template vector_type<double, 2>;
} && requires(
    typename P::template matrix_type<double, 2, 3> a23,
    typename P::template matrix_type<double, 3, 2> a32,
    typename P::template matrix_type<double, 2, 2> sq,
    typename P::template matrix_type<double, 2, 2> sq2,
    typename P::template vector_type<double, 2> v
) {
    { P::multiply(a23, a32) } -> std::same_as<typename P::template matrix_type<double, 2, 2>>;
    { P::multiply(sq, v) } -> std::same_as<typename P::template vector_type<double, 2>>;
    { P::transpose(a23) } -> std::same_as<typename P::template matrix_type<double, 3, 2>>;
    { P::solve(sq, v) } -> std::same_as<typename P::template vector_type<double, 2>>;
    { P::solve(sq, sq2) } -> std::same_as<typename P::template matrix_type<double, 2, 2>>;
    { P::template identity<double, 2>() } -> std::same_as<typename P::template matrix_type<double, 2, 2>>;
    { P::add(a23, a23) } -> std::same_as<typename P::template matrix_type<double, 2, 3>>;
    { P::add(v, v) } -> std::same_as<typename P::template vector_type<double, 2>>;
    { P::subtract(a23, a23) } -> std::same_as<typename P::template matrix_type<double, 2, 3>>;
    { P::subtract(v, v) } -> std::same_as<typename P::template vector_type<double, 2>>;
};

}

#endif
