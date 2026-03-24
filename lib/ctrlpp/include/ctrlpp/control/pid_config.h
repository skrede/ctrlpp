#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_CONFIG_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_CONFIG_H

#include "ctrlpp/types.h"
#include "ctrlpp/pid_policies.h"

#include <tuple>
#include <limits>
#include <cstddef>
#include <type_traits>

namespace ctrlpp {

namespace detail {

// --- Detect whether a policy P has a nested config type ---
// Policies may have either P::config (non-template) or P::template config<Scalar, N> (template).

// Detector for P::template config<Scalar, N>
template <typename P, typename Scalar, std::size_t N, typename = void>
struct has_template_config : std::false_type
{
};

template <typename P, typename Scalar, std::size_t N>
struct has_template_config<P, Scalar, N, std::void_t<typename P::template config<Scalar, N>>> : std::true_type
{
};

// Detector for P::config (non-template)
template <typename P, typename = void>
struct has_plain_config : std::false_type
{
};

template <typename P>
struct has_plain_config<P, std::void_t<typename P::config>> : std::true_type
{
};

// Resolve the config type for a given policy
template <typename P, typename Scalar, std::size_t N, typename = void>
struct policy_config_type;

template <typename P, typename Scalar, std::size_t N>
struct policy_config_type<P, Scalar, N, std::enable_if_t<has_template_config<P, Scalar, N>::value>>
{
    using type = typename P::template config<Scalar, N>;
};

template <typename P, typename Scalar, std::size_t N>
struct policy_config_type<P, Scalar, N, std::enable_if_t<!has_template_config<P, Scalar, N>::value && has_plain_config<P>::value>>
{
    using type = typename P::config;
};

// Check if policy has any config at all
template <typename P, typename Scalar, std::size_t N>
inline constexpr bool has_any_config_v = has_template_config<P, Scalar, N>::value || has_plain_config<P>::value;

// Check if a type appears in a tuple
template <typename T, typename Tuple>
inline constexpr bool tuple_has_v = false;

template <typename T, typename... Ts>
inline constexpr bool tuple_has_v<T, std::tuple<Ts...>> = (std::is_same_v<T, Ts> || ...);

// Build a tuple of config types from policies, filtering out those without config
template <typename Scalar, std::size_t N, typename EnabledList, typename... Policies>
struct policy_configs_builder;

template <typename Scalar, std::size_t N, typename... Collected>
struct policy_configs_builder<Scalar, N, std::tuple<Collected...>>
{
    using type = std::tuple<Collected...>;
};

template <typename Scalar, std::size_t N, typename... Collected, typename P, typename... Rest> requires has_any_config_v<P, Scalar, N>
struct policy_configs_builder<Scalar, N, std::tuple<Collected...>, P, Rest...>
    : policy_configs_builder<Scalar, N, std::tuple<Collected..., typename policy_config_type<P, Scalar, N>::type>, Rest...>
{
};

template <typename Scalar, std::size_t N, typename... Collected, typename P, typename... Rest> requires (!has_any_config_v<P, Scalar, N>)
struct policy_configs_builder<Scalar, N, std::tuple<Collected...>, P, Rest...>
    : policy_configs_builder<Scalar, N, std::tuple<Collected...>, Rest...>
{
};

template <typename Scalar, std::size_t N, typename... Policies>
using policy_configs_tuple_t = typename policy_configs_builder<Scalar, N, std::tuple<>, Policies...>::type;

}

template <typename Scalar, std::size_t NY, typename... Policies>
struct pid_config
{
    using vector_t = Vector<Scalar, NY>;
    using policies_tuple_t = detail::policy_configs_tuple_t<Scalar, NY, Policies...>;

    vector_t kp = vector_t::Zero();
    vector_t ki = vector_t::Zero();
    vector_t kd = vector_t::Zero();
    vector_t output_min = vector_t::Constant(std::numeric_limits<Scalar>::lowest());
    vector_t output_max = vector_t::Constant(std::numeric_limits<Scalar>::max());
    bool derivative_on_error = false;
    vector_t b = vector_t::Constant(Scalar{1});
    vector_t c = vector_t::Constant(Scalar{1});

    policies_tuple_t policies{};

    template <typename P> requires (detail::has_any_config_v<P, Scalar, NY> && detail::tuple_has_v<typename detail::policy_config_type<P, Scalar, NY>::type, policies_tuple_t>)
    constexpr auto &policy()
    {
        using cfg_t = typename detail::policy_config_type<P, Scalar, NY>::type;
        return std::get<cfg_t>(policies);
    }

    template <typename P> requires (detail::has_any_config_v<P, Scalar, NY> && detail::tuple_has_v<typename detail::policy_config_type<P, Scalar, NY>::type, policies_tuple_t>)
    constexpr const auto &policy() const
    {
        using cfg_t = typename detail::policy_config_type<P, Scalar, NY>::type;
        return std::get<cfg_t>(policies);
    }
};

}

#endif
