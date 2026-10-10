#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_GAINS_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_GAINS_H

#include "ctrlpp/control/detail/pid_config_checks.h"

#include "ctrlpp/control/pid_error.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY>
struct pid_gains
{
    Vector<Scalar, NY> kp;
    Vector<Scalar, NY> ki;
    Vector<Scalar, NY> kd;
    Vector<Scalar, NY> kb;
    Vector<Scalar, NY> tf;
};

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

// The default derivative filter time constant is Tf = Td/N, which is kd/(kp*n) in parallel gains.
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3
template <typename Scalar, std::size_t NY, typename... Policies>
auto derivative_filter_times(const pid_config<Scalar, NY, Policies...>& cfg, const pid_gains<Scalar, NY>& gains) -> Vector<Scalar, NY>
{
    const auto& filter = cfg.template policy<deriv_filter>();
    Vector<Scalar, NY> result = Vector<Scalar, NY>::Zero();
    for(std::size_t i = 0; i < NY; ++i)
    {
        const auto row = static_cast<Eigen::Index>(i);
        if(gains.kd[row] != Scalar{0})
            result[row] = filter.tf[i] ? *filter.tf[i] : gains.kd[row] / (gains.kp[row] * filter.n[i]);
    }
    return result;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto compute_pid_gains(const pid_config<Scalar, NY, Policies...>& cfg) -> pid_gains<Scalar, NY>
{
    pid_gains<Scalar, NY> gains{cfg.kp, cfg.ki, cfg.kd, Vector<Scalar, NY>::Zero(), Vector<Scalar, NY>::Zero()};
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
    if constexpr(contains_v<deriv_filter, Policies...>)
        gains.tf = derivative_filter_times(cfg, gains);
    return gains;
}

template <typename Scalar>
auto zero_or_same_sign(Scalar a, Scalar b) -> bool
{
    return a == Scalar{0} || b == Scalar{0} || (a < Scalar{0}) == (b < Scalar{0});
}

// Ti = kp/ki and Td = kd/kp in Tt = sqrt(Ti*Td) are both positive exactly when the nonzero gains share a
// sign, so only then does the default kb exist.
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3
template <typename Scalar, std::size_t NY, typename... Policies>
auto tracking_times_defined(const pid_config<Scalar, NY, Policies...>& cfg, const pid_gains<Scalar, NY>& gains) -> bool
{
    if constexpr(selects_policy<anti_windup, anti_windup<back_calc>, Policies...>())
    {
        const auto& kb = cfg.template policy<anti_windup<back_calc>>().kb;
        for(std::size_t i = 0; i < NY; ++i)
        {
            const auto row = static_cast<Eigen::Index>(i);
            if(kb[i] == Scalar{0} && gains.ki[row] != Scalar{0}
               && !(zero_or_same_sign(gains.ki[row], gains.kp[row]) && zero_or_same_sign(gains.ki[row], gains.kd[row])))
                return false;
        }
    }
    return true;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto derivative_filters_defined(const pid_config<Scalar, NY, Policies...>& cfg, const pid_gains<Scalar, NY>& gains) -> bool
{
    if constexpr(contains_v<deriv_filter, Policies...>)
    {
        const auto& filter = cfg.template policy<deriv_filter>();
        for(std::size_t i = 0; i < NY; ++i)
        {
            const auto row = static_cast<Eigen::Index>(i);
            const Scalar kp = gains.kp[row];
            const Scalar kd = gains.kd[row];
            if(kd != Scalar{0} && !filter.tf[i] && (kp == Scalar{0} || filter.n[i] <= Scalar{0} || !zero_or_same_sign(kd, kp)))
                return false;
        }
    }
    return true;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto derived_gains_representable(const pid_config<Scalar, NY, Policies...>&, const pid_gains<Scalar, NY>& gains) -> bool
{
    if(!gains.kp.allFinite() || !gains.ki.allFinite() || !gains.kd.allFinite() || !gains.kb.allFinite())
        return false;
    if constexpr(contains_v<deriv_filter, Policies...>)
    {
        constexpr Scalar inf = std::numeric_limits<Scalar>::infinity();
        return (gains.kd.array() == Scalar{0} || (gains.tf.array() > Scalar{0} && gains.tf.array() != inf)).all();
    }
    else
        return true;
}

template <typename Scalar, std::size_t NY, typename... Policies>
auto derive_pid_gains(const pid_config<Scalar, NY, Policies...>& cfg) -> ctrlpp::expected<pid_gains<Scalar, NY>, pid_config_error>
{
    if(auto base = check_base_fields(cfg); !base)
        return ctrlpp::unexpected(base.error());
    if(auto policies = check_policy_fields(cfg); !policies)
        return ctrlpp::unexpected(policies.error());
    auto gains = compute_pid_gains(cfg);
    if(!tracking_times_defined(cfg, gains))
        return ctrlpp::unexpected(pid_config_error::undefined_tracking_time);
    if(!derivative_filters_defined(cfg, gains))
        return ctrlpp::unexpected(pid_config_error::undefined_derivative_filter);
    if(!derived_gains_representable(cfg, gains))
        return ctrlpp::unexpected(pid_config_error::unrepresentable_gain);
    return gains;
}

}

#endif
