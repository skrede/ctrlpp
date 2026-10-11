#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_PERFORMANCE_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_PERFORMANCE_H

/// @brief PID performance assessment: IAE, ISE, ITAE, oscillation detection.
///
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3 (IAE, ISE, ITAE)
/// @cite hagglund1995 -- Hagglund, "A control-loop performance monitor", 1995 (oscillation detection)

#include "ctrlpp/types.h"

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <algorithm>

namespace ctrlpp
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_performance_tracker
{
    using vector_t = Vector<Scalar, NY>;
    using config_t = pid_config<Scalar, NY, Policies...>;

public:
    void accumulate(const vector_t& e, Scalar dt, const config_t& cfg)
    {
        if constexpr(detail::has_policy_v<perf_assessment, Policies...>)
        {
            using PA = detail::find_policy_t<perf_assessment, Policies...>;
            m_accumulated_time += dt;

            if constexpr(detail::perf_has_metric_v<IAE, PA>)
                accumulate_iae(e, dt);
            if constexpr(detail::perf_has_metric_v<ISE, PA>)
                accumulate_ise(e, dt);
            if constexpr(detail::perf_has_metric_v<ITAE, PA>)
                accumulate_itae(e, dt);
            if constexpr(detail::perf_has_metric_v<oscillation_detect, PA>)
                accumulate_oscillation(e, dt, cfg);
        }
    }

    template <typename Metric>
    auto metric() const -> const vector_t&
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        static_assert(detail::perf_has_metric_v<Metric, PA>, "Metric type not in perf_assessment pack");

        if constexpr(std::is_same_v<Metric, IAE>)
            return m_iae;
        else if constexpr(std::is_same_v<Metric, ISE>)
            return m_ise;
        else if constexpr(std::is_same_v<Metric, ITAE>)
            return m_itae;
        else if constexpr(std::is_same_v<Metric, oscillation_detect>)
            return m_load_rate;
    }

    auto oscillating(const config_t& cfg) const -> bool
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        static_assert(detail::perf_has_metric_v<oscillation_detect, PA>, "oscillation_detect not in perf_assessment pack");
        return (m_load_rate.array() > cfg.template policy<PA>().load_count_limit).any();
    }

    void reset()
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        if constexpr(detail::perf_has_metric_v<IAE, PA>)
            m_iae = vector_t::Zero();
        if constexpr(detail::perf_has_metric_v<ISE, PA>)
            m_ise = vector_t::Zero();
        if constexpr(detail::perf_has_metric_v<ITAE, PA>)
            m_itae = vector_t::Zero();
        if constexpr(detail::perf_has_metric_v<oscillation_detect, PA>)
        {
            m_crossing_iae = vector_t::Zero();
            m_load_rate = vector_t::Zero();
            m_prev_error_sign = vector_t::Zero();
        }
        m_accumulated_time = Scalar{0};
    }

    auto all_finite() const -> bool
    {
        return m_iae.allFinite()
            && m_ise.allFinite()
            && m_itae.allFinite()
            && m_crossing_iae.allFinite()
            && m_load_rate.allFinite()
            && m_prev_error_sign.allFinite()
            && std::isfinite(m_accumulated_time);
    }

private:
    /// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3 (IAE)
    void accumulate_iae(const vector_t& e, Scalar dt)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
            m_iae[i] += abs_e * dt;
        }
    }

    /// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3 (ISE)
    void accumulate_ise(const vector_t& e, Scalar dt)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            m_ise[i] += e[i] * e[i] * dt;
    }

    /// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3 (ITAE)
    void accumulate_itae(const vector_t& e, Scalar dt)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
            m_itae[i] += m_accumulated_time * abs_e * dt;
        }
    }

    /// @cite hagglund2002 -- Hagglund, "Industrial applications of automatic performance monitoring tools", 2002, Sec. 2.1
    void accumulate_oscillation(const vector_t& e, Scalar dt, const config_t& cfg)
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        const auto& detector = cfg.template policy<PA>();
        for(std::size_t i = 0; i < NY; ++i)
        {
            const auto row = static_cast<Eigen::Index>(i);
            const Scalar sign_e = (e[row] > Scalar{0}) ? Scalar{1} : (e[row] < Scalar{0}) ? Scalar{-1} : Scalar{0};
            const bool crossed = sign_e != Scalar{0} && m_prev_error_sign[row] != Scalar{0} && sign_e != m_prev_error_sign[row];
            const Scalar iae_limit = Scalar{2} * detector.amplitude[i] / detector.ultimate_frequency[i];
            const Scalar load = crossed && m_crossing_iae[row] > iae_limit ? Scalar{1} : Scalar{0};
            const Scalar gamma = forgetting_factor(detector.supervision_periods, detector.ultimate_frequency[i], dt);
            m_crossing_iae[row] = (crossed ? Scalar{0} : m_crossing_iae[row]) + std::abs(e[row]) * dt;
            m_load_rate[row] = gamma * m_load_rate[row] + load;
            if(sign_e != Scalar{0})
                m_prev_error_sign[row] = sign_e;
        }
    }

    // Eq. (2) of hagglund2002 assumes a step much shorter than the window; a step spanning the whole window
    // leaves nothing older, so the factor stops at zero instead of turning negative.
    static auto forgetting_factor(Scalar supervision_periods, Scalar ultimate_frequency, Scalar dt) -> Scalar
    {
        const Scalar supervision_time = supervision_periods * Scalar{2} * std::numbers::pi_v<Scalar> / ultimate_frequency;
        return std::max(Scalar{0}, Scalar{1} - dt / supervision_time);
    }

    vector_t m_iae = vector_t::Zero();
    vector_t m_ise = vector_t::Zero();
    vector_t m_itae = vector_t::Zero();
    vector_t m_crossing_iae = vector_t::Zero();
    vector_t m_load_rate = vector_t::Zero();
    vector_t m_prev_error_sign = vector_t::Zero();
    Scalar m_accumulated_time{0};
};

}

#endif
