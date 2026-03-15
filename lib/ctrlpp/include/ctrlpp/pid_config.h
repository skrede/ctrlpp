#ifndef HPP_GUARD_CPPCTRL_PID_CONFIG_H
#define HPP_GUARD_CPPCTRL_PID_CONFIG_H

#include "ctrlpp/pid_policies.h"
#include "ctrlpp/linalg_policy.h"

#include <cstddef>
#include <limits>
#include <tuple>
#include <type_traits>

namespace ctrlpp {

namespace detail {

// --- Detect whether a policy P has a nested config type ---
// Policies may have either P::config (non-template) or P::template config<Scalar, N> (template).

// Detector for P::template config<Scalar, N>
template<typename P, typename Scalar, std::size_t N, typename = void>
struct has_template_config : std::false_type {};

template<typename P, typename Scalar, std::size_t N>
struct has_template_config<P, Scalar, N,
    std::void_t<typename P::template config<Scalar, N>>> : std::true_type {};

// Detector for P::config (non-template)
template<typename P, typename = void>
struct has_plain_config : std::false_type {};

template<typename P>
struct has_plain_config<P, std::void_t<typename P::config>> : std::true_type {};

// Resolve the config type for a given policy
template<typename P, typename Scalar, std::size_t N, typename = void>
struct policy_config_type;

template<typename P, typename Scalar, std::size_t N>
struct policy_config_type<P, Scalar, N,
    std::enable_if_t<has_template_config<P, Scalar, N>::value>> {
    using type = typename P::template config<Scalar, N>;
};

template<typename P, typename Scalar, std::size_t N>
struct policy_config_type<P, Scalar, N,
    std::enable_if_t<!has_template_config<P, Scalar, N>::value &&
                     has_plain_config<P>::value>> {
    using type = typename P::config;
};

// Check if policy has any config at all
template<typename P, typename Scalar, std::size_t N>
inline constexpr bool has_any_config_v =
    has_template_config<P, Scalar, N>::value || has_plain_config<P>::value;

// Build a tuple of config types from policies, filtering out those without config
template<typename Scalar, std::size_t N, typename EnabledList, typename... Policies>
struct policy_configs_builder;

template<typename Scalar, std::size_t N, typename... Collected>
struct policy_configs_builder<Scalar, N, std::tuple<Collected...>> {
    using type = std::tuple<Collected...>;
};

template<typename Scalar, std::size_t N, typename... Collected, typename P, typename... Rest>
    requires has_any_config_v<P, Scalar, N>
struct policy_configs_builder<Scalar, N, std::tuple<Collected...>, P, Rest...>
    : policy_configs_builder<Scalar, N,
          std::tuple<Collected..., typename policy_config_type<P, Scalar, N>::type>,
          Rest...> {};

template<typename Scalar, std::size_t N, typename... Collected, typename P, typename... Rest>
    requires (!has_any_config_v<P, Scalar, N>)
struct policy_configs_builder<Scalar, N, std::tuple<Collected...>, P, Rest...>
    : policy_configs_builder<Scalar, N, std::tuple<Collected...>, Rest...> {};

template<typename Scalar, std::size_t N, typename... Policies>
using policy_configs_tuple_t =
    typename policy_configs_builder<Scalar, N, std::tuple<>, Policies...>::type;

}

template<LinalgPolicy Policy, typename Scalar, std::size_t NY, typename... Policies>
struct PidConfig {
    using vector_t = typename Policy::template vector_type<Scalar, NY>;
    using policies_tuple_t = detail::policy_configs_tuple_t<Scalar, NY, Policies...>;

    vector_t kp{};
    vector_t ki{};
    vector_t kd{};
    vector_t output_min = []() {
        vector_t v{};
        for (std::size_t i = 0; i < NY; ++i)
            v[i] = std::numeric_limits<Scalar>::lowest();
        return v;
    }();
    vector_t output_max = []() {
        vector_t v{};
        for (std::size_t i = 0; i < NY; ++i)
            v[i] = std::numeric_limits<Scalar>::max();
        return v;
    }();
    bool derivative_on_error = false;
    vector_t b = []() {
        vector_t v{};
        for (std::size_t i = 0; i < NY; ++i)
            v[i] = Scalar{1};
        return v;
    }();
    vector_t c = []() {
        vector_t v{};
        for (std::size_t i = 0; i < NY; ++i)
            v[i] = Scalar{1};
        return v;
    }();

    policies_tuple_t policies{};

    template<typename P>
        requires detail::has_any_config_v<P, Scalar, NY>
    constexpr auto& policy() {
        using cfg_t = typename detail::policy_config_type<P, Scalar, NY>::type;
        return std::get<cfg_t>(policies);
    }

    template<typename P>
        requires detail::has_any_config_v<P, Scalar, NY>
    constexpr const auto& policy() const {
        using cfg_t = typename detail::policy_config_type<P, Scalar, NY>::type;
        return std::get<cfg_t>(policies);
    }
};

}

#endif
