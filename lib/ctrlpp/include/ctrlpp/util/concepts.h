#ifndef HPP_GUARD_CTRLPP_UTIL_CONCEPTS_H
#define HPP_GUARD_CTRLPP_UTIL_CONCEPTS_H

/// @brief Shared C++20 concept constraints used across ctrlpp public APIs.
///
/// Centralises the scalar-type constraint that previously lived as a
/// scattered `static_assert(std::is_floating_point_v<Scalar>, ...)` at each
/// template entry point. Using a named concept at the template head lets the
/// compiler reject malformed instantiations before any body is parsed, gives
/// clearer diagnostics, and declutters class and function bodies.

#include <concepts>

namespace ctrlpp
{

/// @brief Concept satisfied by any floating-point scalar type accepted by
///        the ctrlpp numerical algorithms (currently `float`, `double`,
///        `long double`). Extensions (e.g. extended-precision scalars that
///        model `std::floating_point`) are picked up automatically without
///        further changes.
template <typename S>
concept ctrlpp_floating_scalar = std::floating_point<S>;

}

#endif
