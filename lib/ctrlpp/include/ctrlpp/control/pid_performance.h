#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_PERFORMANCE_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_PERFORMANCE_H

/// @brief PID performance assessment: IAE, ISE, ITAE, oscillation detection.
///
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3

#include "ctrlpp/types.h"

#include "ctrlpp/control/pid_policies.h"

#include <cmath>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_performance_tracker
{
    using vector_t = Vector<Scalar, NY>;

public:
    void accumulate(const vector_t& e, Scalar dt)
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
                accumulate_oscillation(e);
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
            return m_zero_crossings;
    }

    auto oscillating() const -> bool
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        static_assert(detail::perf_has_metric_v<oscillation_detect, PA>, "oscillation_detect not in perf_assessment pack");

        if(m_accumulated_time <= Scalar{0})
            return false;

        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if(m_zero_crossings[i] / m_accumulated_time > static_cast<Scalar>(m_osc_threshold))
                return true;
        }
        return false;
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
            m_zero_crossings = vector_t::Zero();
            m_prev_error_sign = vector_t::Zero();
        }
        m_accumulated_time = Scalar{0};
    }

    void set_first_step(bool first) { m_first_step = first; }

    // Set the zero-crossing rate above which oscillating() reports a limit cycle.
    // The owning controller wires this from perf_assessment::config so the exposed
    // default is the single source of truth.
    void set_oscillation_threshold(Scalar threshold) { m_osc_threshold = threshold; }

    auto all_finite() const -> bool
    {
        return m_iae.allFinite()
            && m_ise.allFinite()
            && m_itae.allFinite()
            && m_zero_crossings.allFinite()
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

    void accumulate_oscillation(const vector_t& e)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            Scalar sign_e = (e[i] > Scalar{0}) ? Scalar{1} : (e[i] < Scalar{0}) ? Scalar{-1} : Scalar{0};
            if(!m_first_step && sign_e != Scalar{0} && m_prev_error_sign[i] != Scalar{0} && sign_e != m_prev_error_sign[i])
                m_zero_crossings[i] += Scalar{1};
            if(sign_e != Scalar{0})
                m_prev_error_sign[i] = sign_e;
        }
    }

    vector_t m_iae = vector_t::Zero();
    vector_t m_ise = vector_t::Zero();
    vector_t m_itae = vector_t::Zero();
    vector_t m_zero_crossings = vector_t::Zero();
    vector_t m_prev_error_sign = vector_t::Zero();
    Scalar m_accumulated_time{0};
    Scalar m_osc_threshold{static_cast<Scalar>(default_crossing_rate_threshold)};
    bool m_first_step{true};
};

}

#endif
