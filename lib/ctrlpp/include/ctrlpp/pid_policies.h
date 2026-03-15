#ifndef HPP_GUARD_CPPCTRL_PID_POLICIES_H
#define HPP_GUARD_CPPCTRL_PID_POLICIES_H

#include "ctrlpp/discretise.h"

#include <array>
#include <cstddef>
#include <type_traits>

namespace ctrlpp {

// --- Anti-windup strategy tags ---

struct BackCalc {};
struct Clamping {};
struct ConditionalIntegration {};

// --- Policy types ---

template<typename Strategy = BackCalc>
struct AntiWindup;

template<>
struct AntiWindup<BackCalc> {
    template<typename Scalar, std::size_t N>
    struct config {
        std::array<Scalar, N> kb{};
    };
};

template<>
struct AntiWindup<Clamping> {
    struct config {};
};

template<>
struct AntiWindup<ConditionalIntegration> {
    template<typename Scalar, std::size_t N>
    struct config {
        std::array<Scalar, N> error_threshold{};
    };
};

struct DerivFilter {
    template<typename Scalar, std::size_t N>
    struct config {
        std::array<Scalar, N> n{};
    };
};

struct SetpointFilter {
    template<typename Scalar, std::size_t N>
    struct config {
        std::array<Scalar, N> tf{};
    };
};

struct PvFilter {
    template<typename Scalar, std::size_t N>
    struct config {
        std::array<Scalar, N> tf{};
    };
};

template<typename Callable = void>
struct FeedForward {
    struct config {
        Callable ff_func;
    };
};

template<>
struct FeedForward<void> {
    struct config {};
};

struct RateLimit {
    template<typename Scalar, std::size_t N>
    struct config {
        std::array<Scalar, N> rate_max{};
    };
};

struct VelocityForm {};
struct IsaForm {};

// --- Performance assessment metric tags ---

struct IAE {};
struct ISE {};
struct ITAE {};

struct OscillationDetect {
    struct config {
        double crossing_rate_threshold{5.0};
    };
};

template<typename... Metrics>
struct PerfAssessment {
    using metrics = std::tuple<Metrics...>;
    struct config {};
};

namespace detail {

// --- contains_v: check if an exact type is in a pack ---

template<typename T, typename... Pack>
struct contains : std::false_type {};

template<typename T, typename First, typename... Rest>
struct contains<T, First, Rest...>
    : std::conditional_t<std::is_same_v<T, First>, std::true_type, contains<T, Rest...>> {};

template<typename T, typename... Pack>
inline constexpr bool contains_v = contains<T, Pack...>::value;

// --- has_policy: check if a template PolicyBase<...> appears in a pack ---

template<template<typename...> class PolicyBase, typename... Policies>
struct has_policy : std::false_type {};

template<template<typename...> class PolicyBase, typename First, typename... Rest>
struct has_policy<PolicyBase, First, Rest...>
    : has_policy<PolicyBase, Rest...> {};

template<template<typename...> class PolicyBase, typename... Args, typename... Rest>
struct has_policy<PolicyBase, PolicyBase<Args...>, Rest...>
    : std::true_type {};

template<template<typename...> class PolicyBase, typename... Policies>
inline constexpr bool has_policy_v = has_policy<PolicyBase, Policies...>::value;

// --- find_policy_t: extract the matching PolicyBase<...> from a pack ---

template<template<typename...> class PolicyBase, typename... Policies>
struct find_policy;

template<template<typename...> class PolicyBase, typename... Args, typename... Rest>
struct find_policy<PolicyBase, PolicyBase<Args...>, Rest...> {
    using type = PolicyBase<Args...>;
};

template<template<typename...> class PolicyBase, typename First, typename... Rest>
struct find_policy<PolicyBase, First, Rest...>
    : find_policy<PolicyBase, Rest...> {};

template<template<typename...> class PolicyBase, typename... Policies>
using find_policy_t = typename find_policy<PolicyBase, Policies...>::type;

// --- perf_has_metric: check if a metric type is in a PerfAssessment's Metrics pack ---

template<typename Metric, typename PA>
struct perf_has_metric : std::false_type {};

template<typename Metric, typename... Metrics>
struct perf_has_metric<Metric, PerfAssessment<Metrics...>>
    : contains<Metric, Metrics...> {};

template<typename Metric, typename PA>
inline constexpr bool perf_has_metric_v = perf_has_metric<Metric, PA>::value;

}

}

#endif
