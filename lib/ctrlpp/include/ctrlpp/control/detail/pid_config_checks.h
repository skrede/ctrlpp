#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_CONFIG_CHECKS_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_CONFIG_CHECKS_H

#include "ctrlpp/control/pid_error.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <optional>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

template <template <typename...> class PolicyBase, typename Selected, typename... Policies>
constexpr auto selects_policy() -> bool
{
    if constexpr(has_policy_v<PolicyBase, Policies...>)
        return std::is_same_v<find_policy_t<PolicyBase, Policies...>, Selected>;
    else
        return false;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto gains_finite(const pid_config<Scalar, NY, Policies...>& cfg) -> bool
{
    if(!cfg.kp.allFinite() || !cfg.kd.allFinite())
        return false;
    // In ISA form ki is the integral time, and +infinity yields the internal gain kp/Ti = 0.
    if constexpr(contains_v<isa_form, Policies...>)
        return !cfg.ki.hasNaN() && (cfg.ki.array() != -std::numeric_limits<Scalar>::infinity()).all();
    else
        return cfg.ki.allFinite();
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto integral_times_nonzero(const pid_config<Scalar, NY, Policies...>& cfg) -> bool
{
    if constexpr(contains_v<isa_form, Policies...>)
        return (cfg.ki.array() != Scalar{0}).all();
    else
        return true;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto output_limits_valid(const pid_config<Scalar, NY, Policies...>& cfg) -> bool
{
    constexpr Scalar inf = std::numeric_limits<Scalar>::infinity();
    const auto lo = cfg.output_min.array();
    const auto hi = cfg.output_max.array();
    return !cfg.output_min.hasNaN() && !cfg.output_max.hasNaN() && (lo <= hi).all() && (lo != inf).all()
        && (hi != -inf).all();
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto check_base_fields(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<void, pid_config_error>
{
    if(!gains_finite(cfg))
        return ctrlpp::unexpected(pid_config_error::non_finite_gain);
    if(!integral_times_nonzero(cfg))
        return ctrlpp::unexpected(pid_config_error::zero_integral_time);
    if(!output_limits_valid(cfg))
        return ctrlpp::unexpected(pid_config_error::invalid_output_limits);
    if(!cfg.b.allFinite() || !cfg.c.allFinite())
        return ctrlpp::unexpected(pid_config_error::non_finite_setpoint_weight);
    return {};
}

template <typename Scalar>
auto filter_time_constant_valid(Scalar tf) -> bool
{
    return tf >= Scalar{0} && tf != std::numeric_limits<Scalar>::infinity();
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto input_filters_valid(const pid_config<Scalar, NY, Policies...>& cfg) -> bool
{
    if constexpr(contains_v<setpoint_filter, Policies...>)
    {
        if(!std::ranges::all_of(cfg.template policy<setpoint_filter>().tf, filter_time_constant_valid<Scalar>))
            return false;
    }
    if constexpr(contains_v<pv_filter, Policies...>)
        return std::ranges::all_of(cfg.template policy<pv_filter>().tf, filter_time_constant_valid<Scalar>);
    else
        return true;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto rate_limit_valid(const pid_config<Scalar, NY, Policies...>& cfg) -> bool
{
    if constexpr(contains_v<rate_limit, Policies...>)
        return std::ranges::all_of(cfg.template policy<rate_limit>().rate_max, [](Scalar r) { return r >= Scalar{0}; });
    else
        return true;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto anti_windup_valid(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<void, pid_config_error>
{
    using back_calc_t = anti_windup<back_calc>;
    using conditional_t = anti_windup<conditional_integration>;
    if constexpr(selects_policy<anti_windup, back_calc_t, Policies...>())
    {
        if(!std::ranges::all_of(cfg.template policy<back_calc_t>().kb, [](Scalar kb) { return std::isfinite(kb); }))
            return ctrlpp::unexpected(pid_config_error::invalid_back_calculation_gain);
    }
    else if constexpr(selects_policy<anti_windup, conditional_t, Policies...>())
    {
        if(std::ranges::any_of(cfg.template policy<conditional_t>().error_threshold, [](Scalar e) { return std::isnan(e); }))
            return ctrlpp::unexpected(pid_config_error::invalid_error_threshold);
    }
    return {};
}

template <typename Scalar>
auto explicit_filter_time_valid(const std::optional<Scalar>& tf) -> bool
{
    return !tf || (*tf > Scalar{0} && std::isfinite(*tf));
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto derivative_filter_valid(const pid_config<Scalar, NY, Policies...>& cfg) -> bool
{
    if constexpr(contains_v<deriv_filter, Policies...>)
    {
        const auto& filter = cfg.template policy<deriv_filter>();
        return std::ranges::all_of(filter.n, [](Scalar n) { return std::isfinite(n); })
            && std::ranges::all_of(filter.tf, explicit_filter_time_valid<Scalar>);
    }
    else
        return true;
}

template <typename Scalar>
auto positive_finite(Scalar value) -> bool
{
    return value > Scalar{0} && std::isfinite(value);
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto oscillation_detector_valid(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<void, pid_config_error>
{
    if constexpr(has_policy_v<perf_assessment, Policies...>)
    {
        using PA = find_policy_t<perf_assessment, Policies...>;
        if constexpr(perf_has_metric_v<oscillation_detect, PA>)
        {
            const auto& detector = cfg.template policy<PA>();
            if(!std::ranges::all_of(detector.amplitude, positive_finite<Scalar>)
               || !std::ranges::all_of(detector.ultimate_frequency, positive_finite<Scalar>))
                return ctrlpp::unexpected(pid_config_error::missing_oscillation_parameters);
            if(!positive_finite(detector.load_count_limit) || !positive_finite(detector.supervision_periods))
                return ctrlpp::unexpected(pid_config_error::invalid_oscillation_supervision);
        }
    }
    return {};
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto check_policy_fields(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<void, pid_config_error>
{
    if(!input_filters_valid(cfg))
        return ctrlpp::unexpected(pid_config_error::invalid_filter_time_constant);
    if(!rate_limit_valid(cfg))
        return ctrlpp::unexpected(pid_config_error::invalid_rate_limit);
    if(auto anti_windup_fields = anti_windup_valid(cfg); !anti_windup_fields)
        return anti_windup_fields;
    if(!derivative_filter_valid(cfg))
        return ctrlpp::unexpected(pid_config_error::invalid_derivative_filter);
    return oscillation_detector_valid(cfg);
}

}

#endif
