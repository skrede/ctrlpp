#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_H

/// @brief Owned C++20 implementation of the ctrlpp expected result type.
///
/// One always-on implementation on every toolchain: ctrlpp does not switch to
/// std::expected even where the standard library ships it. Storage is a raw
/// discriminated union rather than std::variant so the embedded floor stays
/// tight: no <variant>, no bad_variant_access / valueless-by-exception
/// machinery, and both throw sites are gated behind CTRLPP_HAS_EXCEPTIONS so the
/// header compiles clean under -fno-exceptions: value() on an error state falls
/// back to std::abort(), and the rollback that restores a cross-state assignment
/// compiles out together with the exception it exists to catch. Per
/// std::expected semantics operator* and error() are unchecked (precondition on
/// has_value()); only value() is checked.
///
/// A cross-state assignment reinitializes the union, so it carries the same
/// constraints std::expected does -- both members assignable and constructible
/// in the relevant flavor, and at least one of them nothrow-move-constructible.
/// A specialization that cannot meet them has no assignment operator rather than
/// an assignment operator that cannot roll back.
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

#include "ctrlpp/config.h"

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

#if CTRLPP_HAS_EXCEPTIONS
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
#if CTRLPP_HAS_EXCEPTIONS
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
        requires(std::is_copy_assignable_v<T> && std::is_copy_constructible_v<T> && std::is_copy_assignable_v<E> && std::is_copy_constructible_v<E> &&
                 (std::is_nothrow_move_constructible_v<T> || std::is_nothrow_move_constructible_v<E>))
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
        requires(std::is_move_assignable_v<T> && std::is_move_constructible_v<T> && std::is_move_assignable_v<E> && std::is_move_constructible_v<E> &&
                 (std::is_nothrow_move_constructible_v<T> || std::is_nothrow_move_constructible_v<E>))
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
    /// Ends the lifetime of *old_member and constructs *new_member from args in
    /// the same union storage. Precondition: the union currently holds
    /// *old_member, and m_has_value still says so.
    ///
    /// Three branches, so that a construction which throws never leaves the union
    /// holding a member whose lifetime already ended:
    ///  1. the new member cannot throw when built from args -- destroy, then
    ///     construct;
    ///  2. it can throw but moves without throwing -- build a temporary first, so
    ///     a throw leaves the union untouched, then destroy and move in;
    ///  3. neither -- stage the OLD member into a local, destroy it, construct the
    ///     new one, and restore the local into its slot if that throws. The
    ///     assignment operators' constraints guarantee that at least one of the
    ///     two members moves without throwing, which is what makes the staging in
    ///     this branch itself unable to fail.
    ///
    /// Every branch leaves m_has_value alone; the caller flips it only after this
    /// returns, so an unwound assignment cannot leave the flag disagreeing with
    /// the union's real occupant.
    template <typename NewT, typename OldT, typename... Args>
    static constexpr void reinit_member(NewT* new_member, OldT* old_member, Args&&... args)
    {
        if constexpr (std::is_nothrow_constructible_v<NewT, Args...>)
        {
            if constexpr (!std::is_trivially_destructible_v<OldT>)
                std::destroy_at(old_member);
            std::construct_at(new_member, std::forward<Args>(args)...);
        }
        else if constexpr (std::is_nothrow_move_constructible_v<NewT>)
        {
            NewT staged(std::forward<Args>(args)...);
            if constexpr (!std::is_trivially_destructible_v<OldT>)
                std::destroy_at(old_member);
            std::construct_at(new_member, std::move(staged));
        }
        else
        {
            static_assert(std::is_nothrow_move_constructible_v<OldT>, "the assignment constraints guarantee the staged member moves back into the union without throwing");
            OldT staged(std::move(*old_member));
            if constexpr (!std::is_trivially_destructible_v<OldT>)
                std::destroy_at(old_member);
#if CTRLPP_HAS_EXCEPTIONS
            try
            {
                std::construct_at(new_member, std::forward<Args>(args)...);
            }
            catch(...)
            {
                std::construct_at(old_member, std::move(staged));
                throw;
            }
#else
            std::construct_at(new_member, std::forward<Args>(args)...);
#endif
        }
    }

    template <typename Arg>
    constexpr void reinit_as_value(Arg&& arg)
    {
        reinit_member(std::addressof(m_value), std::addressof(m_error), std::forward<Arg>(arg));
        m_has_value = true;
    }

    template <typename Arg>
    constexpr void reinit_as_error(Arg&& arg)
    {
        reinit_member(std::addressof(m_error), std::addressof(m_value), std::forward<Arg>(arg));
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
