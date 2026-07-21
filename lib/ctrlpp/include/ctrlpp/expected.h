#ifndef HPP_GUARD_CTRLPP_EXPECTED_H
#define HPP_GUARD_CTRLPP_EXPECTED_H

/// @brief Vocabulary result type for every fallible ctrlpp operation.
///
/// `ctrlpp::expected` is a single always-on C++20 implementation that ctrlpp
/// owns on every toolchain; it never aliases `std::expected`, so the same
/// storage, accessors, and exception-gated `value()` ship regardless of the
/// standard library vintage. It is a faithful-API result type, not bit-for-bit
/// `std::expected`: the owned special members make it non-trivially-copyable
/// even for trivial `T`/`E`, so triviality is not propagated. Boundary interop
/// with `std::expected` / `std::unexpected` is available through the explicit
/// converting constructors and conversion operators when the standard library
/// provides them.

#include "ctrlpp/detail/expected.h"

#include <utility>
#include <version>
#include <type_traits>

#if defined(__cpp_lib_expected)
    #include <expected>
#endif

namespace ctrlpp
{

// A class template with a deduction guide rather than an alias template: alias-
// template CTAD (P1814) is GCC-only, so `unexpected(err)` through an alias fails
// to compile on Clang. A real class template deduces E on every compiler and
// mirrors the shape of std::unexpected.
template <typename E>
class unexpected
{
    E m_error;

public:
    explicit constexpr unexpected(E err)
        : m_error(std::move(err))
    {
    }

#if defined(__cpp_lib_expected)
    explicit constexpr unexpected(const std::unexpected<E>& other)
        : m_error(other.error())
    {
    }

    explicit constexpr unexpected(std::unexpected<E>&& other)
        : m_error(std::move(other).error())
    {
    }

    explicit constexpr operator std::unexpected<E>() const&
    {
        return std::unexpected<E>(m_error);
    }

    explicit constexpr operator std::unexpected<E>() &&
    {
        return std::unexpected<E>(std::move(m_error));
    }
#endif

    constexpr E& error() & noexcept
    {
        return m_error;
    }
    constexpr const E& error() const& noexcept
    {
        return m_error;
    }
    constexpr E&& error() && noexcept
    {
        return std::move(m_error);
    }
    constexpr const E&& error() const&& noexcept
    {
        return std::move(m_error);
    }
};

template <typename E>
unexpected(E) -> unexpected<E>;

using detail::expected;
using detail::unexpect_t;
using detail::unexpect;

}

#endif
