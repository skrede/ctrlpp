#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_POLICY_TRAITS_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_POLICY_TRAITS_H

#include <type_traits>

namespace ctrlpp::detail
{

template <typename T, typename... Pack>
struct contains : std::false_type
{
};

template <typename T, typename First, typename... Rest>
struct contains<T, First, Rest...> : std::conditional_t<std::is_same_v<T, First>, std::true_type, contains<T, Rest...>>
{
};

template <typename T, typename... Pack>
inline constexpr bool contains_v = contains<T, Pack...>::value;

template <template <typename...> class PolicyBase, typename... Policies>
struct has_policy : std::false_type
{
};

template <template <typename...> class PolicyBase, typename First, typename... Rest>
struct has_policy<PolicyBase, First, Rest...> : has_policy<PolicyBase, Rest...>
{
};

template <template <typename...> class PolicyBase, typename... Args, typename... Rest>
struct has_policy<PolicyBase, PolicyBase<Args...>, Rest...> : std::true_type
{
};

template <template <typename...> class PolicyBase, typename... Policies>
inline constexpr bool has_policy_v = has_policy<PolicyBase, Policies...>::value;

template <template <typename...> class PolicyBase, typename... Policies>
struct find_policy;

template <template <typename...> class PolicyBase, typename... Args, typename... Rest>
struct find_policy<PolicyBase, PolicyBase<Args...>, Rest...>
{
    using type = PolicyBase<Args...>;
};

template <template <typename...> class PolicyBase, typename First, typename... Rest>
struct find_policy<PolicyBase, First, Rest...> : find_policy<PolicyBase, Rest...>
{
};

template <template <typename...> class PolicyBase, typename... Policies>
using find_policy_t = typename find_policy<PolicyBase, Policies...>::type;

}

#endif
