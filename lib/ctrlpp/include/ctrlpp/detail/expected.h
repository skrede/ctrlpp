#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_H

/// @brief Owned C++20 implementation of the ctrlpp expected result type.
///
/// One always-on implementation on every toolchain: ctrlpp does not switch to
/// std::expected even where the standard library ships it. Storage is a raw
/// discriminated union rather than std::variant so the embedded floor stays
/// tight: no <variant>, no bad_variant_access / valueless-by-exception
/// machinery, and the only throw site (value() on an error) is gated behind
/// __cpp_exceptions with a std::abort() fallback so the header compiles clean
/// under -fno-exceptions. Per std::expected semantics operator* and error()
/// are unchecked (precondition on has_value()); only value() is checked.
///
/// Faithful-API, not bit-for-bit std::expected: the owned copy/move special
/// members make this non-trivially-copyable even when T and E are trivial, so
/// triviality is not propagated to the ABI. ctrlpp's T is usually a non-trivial
/// Eigen type, so the cost rarely bites; revisit only if a caller needs a
/// trivially-relocatable expected.
///
/// The unexpected error wrapper and its deduction guide live in the public
/// "ctrlpp/expected.h"; the constructors here take that ctrlpp::unexpected<E>.
/// Explicit converting constructors and conversion operators to and from
/// std::expected are provided under __cpp_lib_expected for boundary interop.
///
/// Monadic operations (and_then, or_else, transform) are intentionally
/// omitted; no caller in the library uses them. Add them here if a future
/// caller needs them.

#include <memory>
#include <cstdlib>
#include <utility>
#include <version>
#include <exception>
#include <type_traits>

#if defined(__cpp_lib_expected)
    #include <expected>
#endif

namespace ctrlpp
{

template <typename E>
class unexpected;

}

namespace ctrlpp::detail
{

struct unexpect_t
{
    explicit unexpect_t() = default;
};

inline constexpr unexpect_t unexpect{};

#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
struct bad_expected_access : std::exception
{
    const char* what() const noexcept override
    {
        return "bad ctrlpp::expected access";
    }
};
#endif

[[noreturn]] inline void on_bad_expected_access()
{
#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
    throw bad_expected_access{};
#else
    std::abort();
#endif
}

template <typename T, typename E>
class expected
{
    bool m_has_value;
    union
    {
        T m_value;
        E m_error;
    };

public:
    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(T val)
        : m_has_value(true)
    {
        std::construct_at(std::addressof(m_value), std::move(val));
    }

    template <typename U, std::enable_if_t<std::is_constructible_v<T, U> && !std::is_same_v<std::remove_cvref_t<U>, expected>, int> = 0>
    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(U&& val)
        : m_has_value(true)
    {
        std::construct_at(std::addressof(m_value), T(std::forward<U>(val)));
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(unexpected<E> err)
        : m_has_value(false)
    {
        std::construct_at(std::addressof(m_error), std::move(err).error());
    }

    template <typename... Args>
    constexpr explicit expected(unexpect_t, Args&&... args)
        : m_has_value(false)
    {
        std::construct_at(std::addressof(m_error), E(std::forward<Args>(args)...));
    }

    constexpr expected(const expected& other)
        : m_has_value(other.m_has_value)
    {
        if (m_has_value)
            std::construct_at(std::addressof(m_value), other.m_value);
        else
            std::construct_at(std::addressof(m_error), other.m_error);
    }

    constexpr expected(expected&& other) noexcept(std::is_nothrow_move_constructible_v<T> && std::is_nothrow_move_constructible_v<E>)
        : m_has_value(other.m_has_value)
    {
        if (m_has_value)
            std::construct_at(std::addressof(m_value), std::move(other.m_value));
        else
            std::construct_at(std::addressof(m_error), std::move(other.m_error));
    }

#if defined(__cpp_lib_expected)
    explicit constexpr expected(const std::expected<T, E>& other)
        : m_has_value(other.has_value())
    {
        if (m_has_value)
            std::construct_at(std::addressof(m_value), *other);
        else
            std::construct_at(std::addressof(m_error), other.error());
    }

    explicit constexpr expected(std::expected<T, E>&& other)
        : m_has_value(other.has_value())
    {
        if (m_has_value)
            std::construct_at(std::addressof(m_value), *std::move(other));
        else
            std::construct_at(std::addressof(m_error), std::move(other).error());
    }

    explicit constexpr operator std::expected<T, E>() const&
    {
        if (m_has_value)
            return std::expected<T, E>(m_value);
        return std::expected<T, E>(std::unexpect, m_error);
    }

    explicit constexpr operator std::expected<T, E>() &&
    {
        if (m_has_value)
            return std::expected<T, E>(std::move(m_value));
        return std::expected<T, E>(std::unexpect, std::move(m_error));
    }
#endif

    constexpr expected& operator=(const expected& other)
    {
        if (m_has_value && other.m_has_value)
            m_value = other.m_value;
        else if (!m_has_value && !other.m_has_value)
            m_error = other.m_error;
        else if (other.m_has_value)
            reinit_as_value(other.m_value);
        else
            reinit_as_error(other.m_error);
        return *this;
    }

    constexpr expected& operator=(expected&& other) noexcept(std::is_nothrow_move_constructible_v<T> && std::is_nothrow_move_assignable_v<T> && std::is_nothrow_move_constructible_v<E> && std::is_nothrow_move_assignable_v<E>)
    {
        if (m_has_value && other.m_has_value)
            m_value = std::move(other.m_value);
        else if (!m_has_value && !other.m_has_value)
            m_error = std::move(other.m_error);
        else if (other.m_has_value)
            reinit_as_value(std::move(other.m_value));
        else
            reinit_as_error(std::move(other.m_error));
        return *this;
    }

    constexpr ~expected()
    {
        if (m_has_value)
        {
            if constexpr (!std::is_trivially_destructible_v<T>)
                std::destroy_at(std::addressof(m_value));
        }
        else if constexpr (!std::is_trivially_destructible_v<E>)
        {
            std::destroy_at(std::addressof(m_error));
        }
    }

    constexpr explicit operator bool() const noexcept
    {
        return m_has_value;
    }
    constexpr bool has_value() const noexcept
    {
        return m_has_value;
    }

    constexpr T& operator*() & noexcept
    {
        return m_value;
    }
    constexpr const T& operator*() const& noexcept
    {
        return m_value;
    }
    constexpr T&& operator*() && noexcept
    {
        return std::move(m_value);
    }
    constexpr const T&& operator*() const&& noexcept
    {
        return std::move(m_value);
    }

    constexpr T* operator->() noexcept
    {
        return std::addressof(m_value);
    }
    constexpr const T* operator->() const noexcept
    {
        return std::addressof(m_value);
    }

    constexpr T& value() &
    {
        if (!m_has_value)
            on_bad_expected_access();
        return m_value;
    }
    constexpr const T& value() const&
    {
        if (!m_has_value)
            on_bad_expected_access();
        return m_value;
    }
    constexpr T&& value() &&
    {
        if (!m_has_value)
            on_bad_expected_access();
        return std::move(m_value);
    }
    constexpr const T&& value() const&&
    {
        if (!m_has_value)
            on_bad_expected_access();
        return std::move(m_value);
    }

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

    template <typename U>
    constexpr T value_or(U&& fallback) const&
    {
        return m_has_value ? m_value : static_cast<T>(std::forward<U>(fallback));
    }

    template <typename U>
    constexpr T value_or(U&& fallback) &&
    {
        return m_has_value ? std::move(m_value) : static_cast<T>(std::forward<U>(fallback));
    }

private:
    template <typename Arg>
    constexpr void reinit_as_value(Arg&& arg)
    {
        if constexpr (!std::is_trivially_destructible_v<E>)
            std::destroy_at(std::addressof(m_error));
        std::construct_at(std::addressof(m_value), std::forward<Arg>(arg));
        m_has_value = true;
    }

    template <typename Arg>
    constexpr void reinit_as_error(Arg&& arg)
    {
        if constexpr (!std::is_trivially_destructible_v<T>)
            std::destroy_at(std::addressof(m_value));
        std::construct_at(std::addressof(m_error), std::forward<Arg>(arg));
        m_has_value = false;
    }
};

template <typename E>
class expected<void, E>
{
    bool m_has_value;
    union
    {
        E m_error;
    };

public:
    constexpr expected() noexcept
        : m_has_value(true)
    {
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(unexpected<E> err)
        : m_has_value(false)
    {
        std::construct_at(std::addressof(m_error), std::move(err).error());
    }

    template <typename... Args>
    constexpr explicit expected(unexpect_t, Args&&... args)
        : m_has_value(false)
    {
        std::construct_at(std::addressof(m_error), E(std::forward<Args>(args)...));
    }

    constexpr expected(const expected& other)
        : m_has_value(other.m_has_value)
    {
        if (!m_has_value)
            std::construct_at(std::addressof(m_error), other.m_error);
    }

    constexpr expected(expected&& other) noexcept(std::is_nothrow_move_constructible_v<E>)
        : m_has_value(other.m_has_value)
    {
        if (!m_has_value)
            std::construct_at(std::addressof(m_error), std::move(other.m_error));
    }

#if defined(__cpp_lib_expected)
    explicit constexpr expected(const std::expected<void, E>& other)
        : m_has_value(other.has_value())
    {
        if (!m_has_value)
            std::construct_at(std::addressof(m_error), other.error());
    }

    explicit constexpr expected(std::expected<void, E>&& other)
        : m_has_value(other.has_value())
    {
        if (!m_has_value)
            std::construct_at(std::addressof(m_error), std::move(other).error());
    }

    explicit constexpr operator std::expected<void, E>() const&
    {
        if (m_has_value)
            return std::expected<void, E>();
        return std::expected<void, E>(std::unexpect, m_error);
    }

    explicit constexpr operator std::expected<void, E>() &&
    {
        if (m_has_value)
            return std::expected<void, E>();
        return std::expected<void, E>(std::unexpect, std::move(m_error));
    }
#endif

    constexpr expected& operator=(const expected& other)
    {
        if (m_has_value == other.m_has_value)
        {
            if (!m_has_value)
                m_error = other.m_error;
        }
        else if (other.m_has_value)
            destroy_error();
        else
            std::construct_at(std::addressof(m_error), other.m_error);
        m_has_value = other.m_has_value;
        return *this;
    }

    constexpr expected& operator=(expected&& other) noexcept(std::is_nothrow_move_constructible_v<E> && std::is_nothrow_move_assignable_v<E>)
    {
        if (m_has_value == other.m_has_value)
        {
            if (!m_has_value)
                m_error = std::move(other.m_error);
        }
        else if (other.m_has_value)
            destroy_error();
        else
            std::construct_at(std::addressof(m_error), std::move(other.m_error));
        m_has_value = other.m_has_value;
        return *this;
    }

    constexpr ~expected()
    {
        destroy_error();
    }

    constexpr explicit operator bool() const noexcept
    {
        return m_has_value;
    }
    constexpr bool has_value() const noexcept
    {
        return m_has_value;
    }

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

private:
    constexpr void destroy_error()
    {
        if (!m_has_value)
            if constexpr (!std::is_trivially_destructible_v<E>)
                std::destroy_at(std::addressof(m_error));
    }
};

}

#endif
