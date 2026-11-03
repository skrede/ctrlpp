#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_VOID_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_VOID_H

#include "ctrlpp/detail/expected_access.h"

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
class expected;

template <typename E>
class expected<void, E>
{
public:
    using value_type = void;
    using error_type = E;
    using unexpected_type = unexpected<E>;

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
        if(!m_has_value)
            std::construct_at(std::addressof(m_error), other.m_error);
    }

    constexpr expected(expected&& other) noexcept(std::is_nothrow_move_constructible_v<E>)
        : m_has_value(other.m_has_value)
    {
        if(!m_has_value)
            std::construct_at(std::addressof(m_error), std::move(other.m_error));
    }

#if defined(__cpp_lib_expected)
    explicit constexpr expected(const std::expected<void, E>& other)
        : m_has_value(other.has_value())
    {
        if(!m_has_value)
            std::construct_at(std::addressof(m_error), other.error());
    }

    explicit constexpr expected(std::expected<void, E>&& other)
        : m_has_value(other.has_value())
    {
        if(!m_has_value)
            std::construct_at(std::addressof(m_error), std::move(other).error());
    }

    explicit constexpr operator std::expected<void, E>() const&
    {
        if(m_has_value)
            return std::expected<void, E>();
        return std::expected<void, E>(std::unexpect, m_error);
    }

    explicit constexpr operator std::expected<void, E>() &&
    {
        if(m_has_value)
            return std::expected<void, E>();
        return std::expected<void, E>(std::unexpect, std::move(m_error));
    }
#endif

    constexpr expected& operator=(const expected& other)
    {
        if(m_has_value == other.m_has_value)
        {
            if(!m_has_value)
                m_error = other.m_error;
        }
        else if(other.m_has_value)
            destroy_error();
        else
            std::construct_at(std::addressof(m_error), other.m_error);
        m_has_value = other.m_has_value;
        return *this;
    }

    constexpr expected& operator=(expected&& other) noexcept(std::is_nothrow_move_constructible_v<E> && std::is_nothrow_move_assignable_v<E>)
    {
        if(m_has_value == other.m_has_value)
        {
            if(!m_has_value)
                m_error = std::move(other.m_error);
        }
        else if(other.m_has_value)
            destroy_error();
        else
            std::construct_at(std::addressof(m_error), std::move(other.m_error));
        m_has_value = other.m_has_value;
        return *this;
    }

    template <typename G>
        requires(std::is_constructible_v<E, const G&> && std::is_assignable_v<E&, const G&>)
    constexpr expected& operator=(const unexpected<G>& err)
    {
        if(m_has_value)
            std::construct_at(std::addressof(m_error), err.error());
        else
            m_error = err.error();
        m_has_value = false;
        return *this;
    }

    template <typename G>
        requires(std::is_constructible_v<E, G> && std::is_assignable_v<E&, G>)
    constexpr expected& operator=(unexpected<G>&& err)
    {
        if(m_has_value)
            std::construct_at(std::addressof(m_error), std::move(err).error());
        else
            m_error = std::move(err).error();
        m_has_value = false;
        return *this;
    }

    constexpr ~expected() { destroy_error(); }

    constexpr explicit operator bool() const noexcept { return m_has_value; }

    constexpr bool has_value() const noexcept { return m_has_value; }

    constexpr E& error() & noexcept { return m_error; }

    constexpr const E& error() const& noexcept { return m_error; }

    constexpr E&& error() && noexcept { return std::move(m_error); }

    constexpr const E&& error() const&& noexcept { return std::move(m_error); }

private:
    bool m_has_value;
    union
    {
        E m_error;
    };

    constexpr void destroy_error()
    {
        if(!m_has_value)
            if constexpr(!std::is_trivially_destructible_v<E>)
                std::destroy_at(std::addressof(m_error));
    }
};

}

#endif
