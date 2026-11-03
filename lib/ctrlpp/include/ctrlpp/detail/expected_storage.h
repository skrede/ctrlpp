#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_STORAGE_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_STORAGE_H

#include "ctrlpp/detail/expected_access.h"

#include "ctrlpp/config.h"

#include <memory>
#include <utility>
#include <version>
#include <type_traits>

#if defined(__cpp_lib_expected)
    #include <expected>
#endif

namespace ctrlpp::detail
{

// A cross-state assignment reinitializes the union, so it carries std::expected's constraints: at least one
// member must move without throwing, or the assignment could not roll back.
template <typename T, typename E, typename GF>
concept error_reinitializable = std::is_constructible_v<E, GF> && std::is_assignable_v<E&, GF> &&
        (std::is_nothrow_constructible_v<E, GF> || std::is_nothrow_move_constructible_v<T> || std::is_nothrow_move_constructible_v<E>);

// A raw union rather than std::variant, so neither <variant> nor its valueless-by-exception machinery is
// pulled in. The user-provided special members make it non-trivially-copyable even for trivial T and E.
template <typename T, typename E>
class expected_storage
{
public:
    template <typename... Args>
    constexpr explicit expected_storage(std::in_place_t, Args&&... args)
        : m_has_value(true)
    {
        std::construct_at(std::addressof(m_value), std::forward<Args>(args)...);
    }

    template <typename... Args>
    constexpr explicit expected_storage(unexpect_t, Args&&... args)
        : m_has_value(false)
    {
        std::construct_at(std::addressof(m_error), std::forward<Args>(args)...);
    }

    constexpr expected_storage(const expected_storage& other)
        : m_has_value(other.m_has_value)
    {
        if(m_has_value)
            std::construct_at(std::addressof(m_value), other.m_value);
        else
            std::construct_at(std::addressof(m_error), other.m_error);
    }

    constexpr expected_storage(expected_storage&& other) noexcept(std::is_nothrow_move_constructible_v<T> && std::is_nothrow_move_constructible_v<E>)
        : m_has_value(other.m_has_value)
    {
        if(m_has_value)
            std::construct_at(std::addressof(m_value), std::move(other.m_value));
        else
            std::construct_at(std::addressof(m_error), std::move(other.m_error));
    }

#if defined(__cpp_lib_expected)
    explicit constexpr expected_storage(const std::expected<T, E>& other)
        : m_has_value(other.has_value())
    {
        if(m_has_value)
            std::construct_at(std::addressof(m_value), *other);
        else
            std::construct_at(std::addressof(m_error), other.error());
    }

    explicit constexpr expected_storage(std::expected<T, E>&& other)
        : m_has_value(other.has_value())
    {
        if(m_has_value)
            std::construct_at(std::addressof(m_value), *std::move(other));
        else
            std::construct_at(std::addressof(m_error), std::move(other).error());
    }
#endif

    constexpr expected_storage& operator=(const expected_storage& other)
        requires(std::is_copy_assignable_v<T> && std::is_copy_constructible_v<T> && std::is_copy_assignable_v<E> && std::is_copy_constructible_v<E> &&
                (std::is_nothrow_move_constructible_v<T> || std::is_nothrow_move_constructible_v<E>))
    {
        if(m_has_value && other.m_has_value)
            m_value = other.m_value;
        else if(!m_has_value && !other.m_has_value)
            m_error = other.m_error;
        else if(other.m_has_value)
            reinit_as_value(other.m_value);
        else
            reinit_as_error(other.m_error);
        return *this;
    }

    constexpr expected_storage& operator=(expected_storage&& other) noexcept(
            std::is_nothrow_move_constructible_v<T> && std::is_nothrow_move_assignable_v<T> && std::is_nothrow_move_constructible_v<E> && std::is_nothrow_move_assignable_v<E>)
        requires(std::is_move_assignable_v<T> && std::is_move_constructible_v<T> && std::is_move_assignable_v<E> && std::is_move_constructible_v<E> &&
                (std::is_nothrow_move_constructible_v<T> || std::is_nothrow_move_constructible_v<E>))
    {
        if(m_has_value && other.m_has_value)
            m_value = std::move(other.m_value);
        else if(!m_has_value && !other.m_has_value)
            m_error = std::move(other.m_error);
        else if(other.m_has_value)
            reinit_as_value(std::move(other.m_value));
        else
            reinit_as_error(std::move(other.m_error));
        return *this;
    }

    constexpr ~expected_storage()
    {
        if(m_has_value)
        {
            if constexpr(!std::is_trivially_destructible_v<T>)
                std::destroy_at(std::addressof(m_value));
        }
        else if constexpr(!std::is_trivially_destructible_v<E>)
        {
            std::destroy_at(std::addressof(m_error));
        }
    }

protected:
    bool m_has_value;
    union
    {
        T m_value;
        E m_error;
    };

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

private:
    // Replaces the union's occupant without ever leaving it holding a member whose lifetime has ended, should
    // the construction throw. m_has_value is left alone; the caller flips it only after this returns.
    template <typename NewT, typename OldT, typename... Args>
    static constexpr void reinit_member(NewT* new_member, OldT* old_member, Args&&... args)
    {
        if constexpr(std::is_nothrow_constructible_v<NewT, Args...>)
            replace_member(new_member, old_member, std::forward<Args>(args)...);
        else if constexpr(std::is_nothrow_move_constructible_v<NewT>)
        {
            NewT staged(std::forward<Args>(args)...);
            replace_member(new_member, old_member, std::move(staged));
        }
        else
            replace_member_restoring(new_member, old_member, std::forward<Args>(args)...);
    }

    template <typename NewT, typename OldT, typename... Args>
    static constexpr void replace_member(NewT* new_member, OldT* old_member, Args&&... args)
    {
        if constexpr(!std::is_trivially_destructible_v<OldT>)
            std::destroy_at(old_member);
        std::construct_at(new_member, std::forward<Args>(args)...);
    }

    template <typename NewT, typename OldT, typename... Args>
    static constexpr void replace_member_restoring(NewT* new_member, OldT* old_member, Args&&... args)
    {
        static_assert(std::is_nothrow_move_constructible_v<OldT>, "the assignment constraints guarantee the staged member moves back into the union without throwing");
#if CTRLPP_HAS_EXCEPTIONS
        OldT staged(std::move(*old_member));
        if constexpr(!std::is_trivially_destructible_v<OldT>)
            std::destroy_at(old_member);
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
        replace_member(new_member, old_member, std::forward<Args>(args)...);
#endif
    }
};

}

#endif
