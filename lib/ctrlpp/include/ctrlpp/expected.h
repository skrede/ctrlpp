#ifndef HPP_GUARD_CTRLPP_EXPECTED_H
#define HPP_GUARD_CTRLPP_EXPECTED_H

#include "ctrlpp/detail/expected.h"

#include <utility>
#include <version>
#include <type_traits>

#if defined(__cpp_lib_expected)
    #include <expected>
#endif

namespace ctrlpp
{

// A class template with a deduction guide, not an alias template: only GCC implements alias-template CTAD
// (P1814), so `unexpected(err)` through an alias fails to compile on Clang.
template <typename E>
class unexpected
{
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

    explicit constexpr operator std::unexpected<E>() const& { return std::unexpected<E>(m_error); }

    explicit constexpr operator std::unexpected<E>() && { return std::unexpected<E>(std::move(m_error)); }
#endif

    constexpr E& error() & noexcept { return m_error; }

    constexpr const E& error() const& noexcept { return m_error; }

    constexpr E&& error() && noexcept { return std::move(m_error); }

    constexpr const E&& error() const&& noexcept { return std::move(m_error); }

private:
    E m_error;
};

template <typename E>
unexpected(E) -> unexpected<E>;

using detail::expected;
using detail::unexpect_t;
using detail::unexpect;

}

#endif
