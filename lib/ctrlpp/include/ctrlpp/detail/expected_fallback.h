#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_FALLBACK_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_FALLBACK_H

/// @brief Hand-rolled C++20 fallback for std::expected.
///
/// Signature-equal to the std::expected subset ctrlpp uses: construction from
/// T, from a convertible U, from unexpected<E>, and in-place unexpect_t tag
/// construction; has_value(), explicit operator bool, ref-qualified operator*,
/// operator->, value(), error(), and value_or(). The expected<void, E>
/// specialization covers fallible operations without a payload.
///
/// Storage is a raw discriminated union rather than std::variant so the
/// embedded floor stays tight: no <variant>, no bad_variant_access /
/// valueless-by-exception machinery, and the only throw site (value() on an
/// error) is gated behind __cpp_exceptions with a std::abort() fallback so the
/// header compiles clean under -fno-exceptions. Per std::expected semantics
/// operator* and error() are unchecked (precondition on has_value()); only
/// value() is checked.
///
/// Monadic operations (and_then, or_else, transform) are intentionally
/// omitted; no caller in the library uses them. Add them here if a future
/// caller needs them.

#include <memory>
#include <cstdlib>
#include <utility>
#include <exception>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename E>
struct unexpected
{
    E value;

    explicit constexpr unexpected(E e)
        : value(std::move(e))
    {
    }
};

struct unexpect_t
{
    explicit unexpect_t() = default;
};

inline constexpr unexpect_t unexpect{};

struct bad_expected_access : std::exception
{
    const char* what() const noexcept override
    {
        return "bad ctrlpp::expected access";
    }
};

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
        std::construct_at(std::addressof(m_error), std::move(err.value));
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
        std::construct_at(std::addressof(m_error), std::move(err.value));
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
