#ifndef HPP_GUARD_CTRLPP_EXPECTED_H
#define HPP_GUARD_CTRLPP_EXPECTED_H

/// @brief Vocabulary result type for every fallible ctrlpp operation.
///
/// `ctrlpp::expected` resolves to `std::expected` when the standard library
/// provides it (C++23 and later) and to an in-library C++20 fallback with an
/// identical call surface otherwise. Callers cannot tell which target is
/// active; both spell construction, `has_value()`, `operator*`, `value()`,
/// and `error()` identically.

#include "ctrlpp/detail/expected.h"

#include <utility>
#include <type_traits>

namespace ctrlpp
{

template <typename T, typename E>
using expected = detail::expected<T, E>;

// A function rather than an alias template: alias-template CTAD (P1814) lets
// GCC deduce E from `unexpected(err)`, but Clang does not implement it, so the
// alias form fails to compile there. The function deduces E on every compiler.
template <typename E>
[[nodiscard]] constexpr detail::unexpected<std::remove_cvref_t<E>> unexpected(E&& error)
{
    return detail::unexpected<std::remove_cvref_t<E>>(std::forward<E>(error));
}

using detail::unexpect_t;
using detail::unexpect;

}

#endif
