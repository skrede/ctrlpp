#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_GAINS_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_GAINS_H

#include "ctrlpp/control/pid_error.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY>
struct pid_gains
{
    Vector<Scalar, NY> kp;
    Vector<Scalar, NY> ki;
    Vector<Scalar, NY> kd;
    Vector<Scalar, NY> kb;
};

template <template <typename...> class PolicyBase, typename Selected, typename... Policies>
constexpr auto selects_policy() -> bool
{
    if constexpr(has_policy_v<PolicyBase, Policies...>)
        return std::is_same_v<find_policy_t<PolicyBase, Policies...>, Selected>;
    else
        return false;
}

// kb multiplies the saturation error (output units) into the integrator rate (output per time), so it is
// 1/Tt with the tracking time Tt = sqrt(Ti*Td) = sqrt(kd/ki) in parallel gains, or Ti = kp/ki without
// derivative action. A pure-I channel has no reference for a tracking time, so back-calculation is off.
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3
template <typename Scalar, std::size_t NY>
auto default_back_calc_gain(const pid_gains<Scalar, NY>& gains, Eigen::Index i) -> Scalar
{
    if(gains.kd[i] != Scalar{0})
        return std::sqrt(gains.ki[i] / gains.kd[i]);
    if(gains.kp[i] != Scalar{0})
        return gains.ki[i] / gains.kp[i];
    return Scalar{0};
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto back_calc_gains(const pid_config<Scalar, NY, Policies...>& cfg, const pid_gains<Scalar, NY>& gains) -> Vector<Scalar, NY>
{
    const auto& kb = cfg.template policy<anti_windup<back_calc>>().kb;
    Vector<Scalar, NY> result;
    for(std::size_t i = 0; i < NY; ++i)
    {
        const auto row = static_cast<Eigen::Index>(i);
        result[row] = kb[i] != Scalar{0} ? kb[i] : default_back_calc_gain(gains, row);
    }
    return result;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto compute_pid_gains(const pid_config<Scalar, NY, Policies...>& cfg) -> pid_gains<Scalar, NY>
{
    pid_gains<Scalar, NY> gains{cfg.kp, cfg.ki, cfg.kd, Vector<Scalar, NY>::Zero()};
    if constexpr(contains_v<isa_form, Policies...>)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            gains.ki[i] = cfg.ki[i] != Scalar{0} ? cfg.kp[i] / cfg.ki[i] : Scalar{0};
            gains.kd[i] = cfg.kp[i] * cfg.kd[i];
        }
    }
    if constexpr(selects_policy<anti_windup, anti_windup<back_calc>, Policies...>())
        gains.kb = back_calc_gains(cfg, gains);
    return gains;
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

template <typename Scalar, std::size_t NY, typename... Policies>
auto check_policy_fields(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<void, pid_config_error>
{
    if(!input_filters_valid(cfg))
        return ctrlpp::unexpected(pid_config_error::invalid_filter_time_constant);
    if(!rate_limit_valid(cfg))
        return ctrlpp::unexpected(pid_config_error::invalid_rate_limit);
    if(auto anti_windup_fields = anti_windup_valid(cfg); !anti_windup_fields)
        return anti_windup_fields;
    if constexpr(contains_v<deriv_filter, Policies...>)
    {
        if(!std::ranges::all_of(cfg.template policy<deriv_filter>().n, [](Scalar n) { return std::isfinite(n); }))
            return ctrlpp::unexpected(pid_config_error::invalid_derivative_filter);
    }
    return {};
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto derive_pid_gains(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<pid_gains<Scalar, NY>, pid_config_error>
{
    if(auto base = check_base_fields(cfg); !base)
        return ctrlpp::unexpected(base.error());
    if(auto policies = check_policy_fields(cfg); !policies)
        return ctrlpp::unexpected(policies.error());
    auto gains = compute_pid_gains(cfg);
    if(!gains.kp.allFinite() || !gains.ki.allFinite() || !gains.kd.allFinite())
        return ctrlpp::unexpected(pid_config_error::unrepresentable_gain);
    return gains;
}

}

#endif
