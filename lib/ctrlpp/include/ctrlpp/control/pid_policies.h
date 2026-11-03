#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_POLICIES_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_POLICIES_H

// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3-4

#include "ctrlpp/control/detail/pid_policy_traits.h"

#include "ctrlpp/model/discretize.h"

#include <array>
#include <tuple>
#include <cstddef>
#include <optional>
#include <type_traits>

namespace ctrlpp
{

struct back_calc
{
};

struct clamping
{
};

struct conditional_integration
{
};

template <typename Strategy = back_calc>
struct anti_windup;

template <>
struct anti_windup<back_calc>
{
    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> kb{};
    };
};

template <>
struct anti_windup<clamping>
{
    struct config
    {
    };
};

template <>
struct anti_windup<conditional_integration>
{
    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> error_threshold{};
    };
};

struct deriv_filter
{
    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> n{};
        std::array<std::optional<Scalar>, N> tf{};
    };
};

struct setpoint_filter
{
    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> tf{};
    };
};

struct pv_filter
{
    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> tf{};
    };
};

template <typename Callable = void>
struct feed_forward
{
    struct config
    {
        Callable ff_func;
    };
};

template <>
struct feed_forward<void>
{
    struct config
    {
    };
};

struct rate_limit
{
    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> rate_max{};
    };
};

struct velocity_form
{
};

struct isa_form
{
};

struct IAE
{
};

struct ISE
{
};

struct ITAE
{
};

struct oscillation_detect
{
};

template <typename... Metrics>
struct perf_assessment
{
    using metrics = std::tuple<Metrics...>;

    template <typename Scalar, std::size_t N>
    struct config
    {
        std::array<Scalar, N> amplitude{};          // error units; required by oscillation_detect
        std::array<Scalar, N> ultimate_frequency{}; // rad/s; required by oscillation_detect
        // @cite hagglund1995 -- Hagglund, "A control-loop performance monitor", 1995 (n_lim, via hagglund2002 Sec. 2.1)
        Scalar load_count_limit{10};
        // @cite hagglund1995 -- T_sup = 5 n_lim T_u, in ultimate periods (via hagglund2002 Sec. 2.1)
        Scalar supervision_periods{50};
    };
};

namespace detail
{

template <typename Metric, typename PA>
struct perf_has_metric : std::false_type
{
};

template <typename Metric, typename... Metrics>
struct perf_has_metric<Metric, perf_assessment<Metrics...>> : contains<Metric, Metrics...>
{
};

template <typename Metric, typename PA>
inline constexpr bool perf_has_metric_v = perf_has_metric<Metric, PA>::value;

}

}

#endif
