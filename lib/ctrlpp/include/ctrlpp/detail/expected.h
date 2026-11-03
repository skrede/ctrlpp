#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_H

#include "ctrlpp/detail/expected_void.h"
#include "ctrlpp/detail/expected_access.h"
#include "ctrlpp/detail/expected_storage.h"

#include <memory>
#include <utility>
#include <version>
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

template <typename T, typename E>
class expected : private expected_storage<T, E>
{
    using storage = expected_storage<T, E>;

public:
    using value_type = T;
    using error_type = E;
    using unexpected_type = unexpected<E>;

    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(T val)
        : storage(std::in_place, std::move(val))
    {
    }

    template <typename U, std::enable_if_t<std::is_constructible_v<T, U> && !std::is_same_v<std::remove_cvref_t<U>, expected>, int> = 0>
    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(U&& val)
        : storage(std::in_place, T(std::forward<U>(val)))
    {
    }

    // NOLINTNEXTLINE(google-explicit-constructor)
    constexpr expected(unexpected<E> err)
        : storage(unexpect, std::move(err).error())
    {
    }

    template <typename... Args>
        requires std::is_constructible_v<T, Args...>
    constexpr explicit expected(std::in_place_t, Args&&... args)
        : storage(std::in_place, std::forward<Args>(args)...)
    {
    }

    template <typename... Args>
    constexpr explicit expected(unexpect_t, Args&&... args)
        : storage(unexpect, E(std::forward<Args>(args)...))
    {
    }

#if defined(__cpp_lib_expected)
    explicit constexpr expected(const std::expected<T, E>& other)
        : storage(other)
    {
    }

    explicit constexpr expected(std::expected<T, E>&& other)
        : storage(std::move(other))
    {
    }

    explicit constexpr operator std::expected<T, E>() const&
    {
        if(m_has_value)
            return std::expected<T, E>(m_value);
        return std::expected<T, E>(std::unexpect, m_error);
    }

    explicit constexpr operator std::expected<T, E>() &&
    {
        if(m_has_value)
            return std::expected<T, E>(std::move(m_value));
        return std::expected<T, E>(std::unexpect, std::move(m_error));
    }
#endif

    // [expected.object.assign]: the error is built in the union itself, so a factory that constructs in place and then fails holds no second T.
    template <typename G>
        requires error_reinitializable<T, E, const G&>
    constexpr expected& operator=(const unexpected<G>& err)
    {
        if(m_has_value)
            reinit_as_error(err.error());
        else
            m_error = err.error();
        return *this;
    }

    template <typename G>
        requires error_reinitializable<T, E, G>
    constexpr expected& operator=(unexpected<G>&& err)
    {
        if(m_has_value)
            reinit_as_error(std::move(err).error());
        else
            m_error = std::move(err).error();
        return *this;
    }

    constexpr explicit operator bool() const noexcept { return m_has_value; }

    constexpr bool has_value() const noexcept { return m_has_value; }

    constexpr T& operator*() & noexcept { return m_value; }

    constexpr const T& operator*() const& noexcept { return m_value; }

    constexpr T&& operator*() && noexcept { return std::move(m_value); }

    constexpr const T&& operator*() const&& noexcept { return std::move(m_value); }

    constexpr T* operator->() noexcept { return std::addressof(m_value); }

    constexpr const T* operator->() const noexcept { return std::addressof(m_value); }

    // Nothing in the library calls value(); library code branches on has_value(), so no library path reaches
    // the throw, or the abort that replaces it, that value() is required to have.
    constexpr T& value() &
    {
        if(!m_has_value)
            on_bad_expected_access();
        return m_value;
    }

    constexpr const T& value() const&
    {
        if(!m_has_value)
            on_bad_expected_access();
        return m_value;
    }

    constexpr T&& value() &&
    {
        if(!m_has_value)
            on_bad_expected_access();
        return std::move(m_value);
    }

    constexpr const T&& value() const&&
    {
        if(!m_has_value)
            on_bad_expected_access();
        return std::move(m_value);
    }

    constexpr E& error() & noexcept { return m_error; }

    constexpr const E& error() const& noexcept { return m_error; }

    constexpr E&& error() && noexcept { return std::move(m_error); }

    constexpr const E&& error() const&& noexcept { return std::move(m_error); }

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
    using storage::m_error;
    using storage::m_value;
    using storage::m_has_value;
    using storage::reinit_as_error;
};

}

#endif
