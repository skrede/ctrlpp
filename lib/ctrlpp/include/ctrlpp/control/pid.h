#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

/// @brief Policy-based PID controller with compile-time feature composition.
///
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/types.h"

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"
#include "ctrlpp/control/pid_performance.h"

#include <cmath>
#include <limits>

namespace ctrlpp
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid
{
public:
    using config_type = pid_config<Scalar, NY, Policies...>;
    using vector_t = Vector<Scalar, NY>;

    explicit pid(const config_type& cfg) : m_cfg{cfg}
    {
        compute_internal_gains(cfg);
        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
            initialize_back_calc_gains();
        initialize_perf_config();
    }

    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt) -> vector_t
    {
        if(dt <= Scalar{0})
            return m_prev_output;

        auto filtered_sp = apply_setpoint_filter(sp, dt);
        auto filtered_meas = apply_pv_filter(meas, dt);
        auto e = (filtered_sp - filtered_meas).eval();
        m_perf.accumulate(e, dt);

        if constexpr(detail::contains_v<velocity_form, Policies...>)
            return compute_velocity_form(e, sp, filtered_sp, filtered_meas, dt);
        else
            return compute_position_form(e, sp, filtered_sp, filtered_meas, dt);
    }

    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt, const vector_t& tracking_signal) -> vector_t
    {
        auto u = compute(sp, meas, dt);
        if(dt <= Scalar{0})
            return u;
        if constexpr(!detail::contains_v<velocity_form, Policies...>)
        {
            auto non_integral = (u - m_integral).eval();
            m_integral = (tracking_signal - non_integral).eval();
        }
        return u;
    }

    void set_params(const config_type& new_cfg)
    {
        // The integral state is stored in output units (each increment is ki*e*dt and
        // enters the output directly), so its contribution to the output is already
        // continuous across a gain change. Leaving the state untouched gives true
        // bumpless transfer; rescaling it would inject a bump instead of removing one.
        m_cfg = new_cfg;
        compute_internal_gains(new_cfg);
        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
            initialize_back_calc_gains();
        initialize_perf_config();
    }

    const vector_t& error() const { return m_prev_error; }
    const vector_t& integral() const { return m_integral; }
    const config_type& params() const { return m_cfg; }
    bool saturated() const { return m_saturated; }

    void reset()
    {
        m_integral = vector_t::Zero();
        m_prev_error = vector_t::Zero();
        m_prev_meas = vector_t::Zero();
        m_prev_sp = vector_t::Zero();
        m_prev_output = vector_t::Zero();
        m_accumulated_output = vector_t::Zero();
        m_prev_ff = vector_t::Zero();
        m_prev_prev_error = vector_t::Zero();
        if constexpr(detail::contains_v<setpoint_filter, Policies...>)
            m_filtered_sp = vector_t::Zero();
        if constexpr(detail::contains_v<pv_filter, Policies...>)
            m_filtered_meas = vector_t::Zero();
        if constexpr(detail::contains_v<deriv_filter, Policies...>)
            m_prev_deriv_filtered = vector_t::Zero();
        m_first_step = true;
        m_integral_frozen = false;
        m_saturated = false;
        if constexpr(detail::has_policy_v<perf_assessment, Policies...>)
        {
            m_perf.reset();
            m_perf.set_first_step(true);
        }
    }

    void freeze_integral(bool freeze = true) { m_integral_frozen = freeze; }
    void set_integral(const vector_t& val) { m_integral = val; }

    template <typename Metric>
    auto metric() const -> const vector_t&
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_perf.template metric<Metric>();
    }

    auto oscillating() const -> bool
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_perf.oscillating();
    }

    void reset_metrics()
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        m_perf.reset();
    }

private:
    void compute_internal_gains(const config_type& cfg)
    {
        if constexpr(detail::contains_v<isa_form, Policies...>)
        {
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(cfg.ki[i] != Scalar{0})
                    m_ki[i] = cfg.kp[i] / cfg.ki[i];
                else
                    m_ki[i] = Scalar{0};
                m_kd[i] = cfg.kp[i] * cfg.kd[i];
                m_kp[i] = cfg.kp[i];
            }
        }
        else
        {
            m_kp = cfg.kp;
            m_ki = cfg.ki;
            m_kd = cfg.kd;
        }
    }

    void initialize_back_calc_gains()
    {
        using AW = detail::find_policy_t<anti_windup, Policies...>;
        if constexpr(std::is_same_v<AW, anti_windup<back_calc>>)
        {
            const auto& aw_cfg = m_cfg.template policy<AW>();
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(aw_cfg.kb[static_cast<std::size_t>(i)] != Scalar{0})
                    kb_[i] = aw_cfg.kb[static_cast<std::size_t>(i)];
                else
                {
                    // The back-calculation gain multiplies the saturation error
                    // (u_sat - u_raw, in output units) into the integrator whose rate
                    // is in output-per-time, so kb must carry units of 1/time. The
                    // Astrom tracking-time-constant default sets kb = 1/Tt with
                    // Tt = sqrt(Ti*Td); in the internal parallel gains Ti = kp/ki and
                    // Td = kd/kp, so 1/Tt = sqrt(ki/kd). With no derivative action the
                    // tracking time collapses to Ti, giving the fallback kb = ki/kp.
                    if(m_kd[i] != Scalar{0})
                        kb_[i] = std::sqrt(m_ki[i] / m_kd[i]);
                    else if(m_kp[i] != Scalar{0})
                        kb_[i] = m_ki[i] / m_kp[i];
                    else
                        // Pure-I controller: no proportional or derivative reference for
                        // a tracking time, so disable back-calculation (kb = 0) rather
                        // than divide by zero.
                        kb_[i] = Scalar{0};
                }
            }
        }
    }

    void initialize_perf_config()
    {
        if constexpr(detail::has_policy_v<perf_assessment, Policies...>)
        {
            using PA = detail::find_policy_t<perf_assessment, Policies...>;
            m_perf.set_oscillation_threshold(
                static_cast<Scalar>(m_cfg.template policy<PA>().crossing_rate_threshold));
        }
    }

    auto apply_setpoint_filter(const vector_t& sp, Scalar dt) -> vector_t
    {
        if constexpr(detail::contains_v<setpoint_filter, Policies...>)
        {
            const auto& tf_sp = m_cfg.template policy<setpoint_filter>().tf;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(tf_sp[static_cast<std::size_t>(i)] > Scalar{0})
                {
                    auto alpha = tf_sp[static_cast<std::size_t>(i)] / (tf_sp[static_cast<std::size_t>(i)] + dt);
                    m_filtered_sp[i] = alpha * m_filtered_sp[i] + (Scalar{1} - alpha) * sp[i];
                }
                else
                    m_filtered_sp[i] = sp[i];
            }
            return m_filtered_sp;
        }
        else
            return sp;
    }

    auto apply_pv_filter(const vector_t& meas, Scalar dt) -> vector_t
    {
        if constexpr(detail::contains_v<pv_filter, Policies...>)
        {
            const auto& tf_pv = m_cfg.template policy<pv_filter>().tf;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(tf_pv[static_cast<std::size_t>(i)] > Scalar{0})
                {
                    auto alpha = tf_pv[static_cast<std::size_t>(i)] / (tf_pv[static_cast<std::size_t>(i)] + dt);
                    m_filtered_meas[i] = alpha * m_filtered_meas[i] + (Scalar{1} - alpha) * meas[i];
                }
                else
                    m_filtered_meas[i] = meas[i];
            }
            return m_filtered_meas;
        }
        else
            return meas;
    }

    auto compute_velocity_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        auto dp = m_kp.cwiseProduct(e - m_prev_error).eval();
        auto di = (m_ki.cwiseProduct(e) * dt).eval();
        auto d_num = (e - m_prev_error * Scalar{2} + m_prev_prev_error).eval();
        auto dd = m_kd.cwiseProduct(d_num / dt).eval();
        auto delta_u = (dp + di + dd).eval();
        delta_u = apply_feed_forward_velocity(delta_u, sp, dt);
        // Clamp the accumulated output to the output limits, not the raw increment,
        // then emit the increment that reaches the clamped accumulated output. Clamping
        // the increment itself would forbid motion against an asymmetric limit (for
        // example output_min = 0 would make the output monotone non-decreasing);
        // clamping the accumulated output keeps it in range while still letting it rise
        // and fall.
        auto target = (m_accumulated_output + delta_u).eval();
        auto clamped = target.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        auto emitted = (clamped - m_accumulated_output).eval();
        m_accumulated_output = clamped;
        update_state(e, filtered_meas, filtered_sp, emitted);
        return emitted;
    }

    auto apply_feed_forward_velocity(vector_t delta_u, const vector_t& sp, Scalar dt) -> vector_t
    {
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                auto ff = m_cfg.template policy<ff_policy_t>().ff_func(sp, dt).eval();
                // The velocity form emits increments, so inject the change in the
                // feed-forward level, not its absolute value; a constant feed-forward
                // then contributes nothing to the increment and the actuator does not
                // drift.
                auto delta_ff = (ff - m_prev_ff).eval();
                m_prev_ff = ff;
                return (delta_u + delta_ff).eval();
            }
        }
        return delta_u;
    }

    /// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3-4
    auto compute_position_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        auto p = compute_proportional_term(filtered_sp, filtered_meas);
        auto [integral_increment, updated_integral] = compute_integral_term(e, dt);
        auto d = compute_derivative_term(filtered_sp, filtered_meas, dt);
        // Unconstrained control command, before rate limiting and output saturation.
        auto u_unconstrained = compute_raw_output(p, updated_integral, d, sp, dt);
        auto u_limited = apply_rate_limit(u_unconstrained, dt);
        auto u_sat = u_limited.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        // saturated() reports output-limit saturation, the clamp against
        // output_min/output_max.
        m_saturated = (u_sat.array() != u_limited.array()).any();
        // Anti-windup feeds back against the fully unconstrained command, so both the
        // rate limiter and the output saturation contribute. Feeding back only the
        // post-rate-limit value would let the integrator wind up freely along a
        // rate-limited ramp, where the applied output already equals the rate-limited
        // command.
        apply_anti_windup(u_sat, u_unconstrained, e, integral_increment, dt);
        update_state(e, filtered_meas, filtered_sp, u_sat);
        return u_sat;
    }

    /// @cite astrom2006 -- Ch. 3.5 (setpoint weighting b parameter)
    auto compute_proportional_term(const vector_t& filtered_sp, const vector_t& filtered_meas) -> vector_t
    {
        vector_t ep;
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            ep[i] = m_cfg.b[i] * filtered_sp[i] - filtered_meas[i];
        return m_kp.cwiseProduct(ep).eval();
    }

    /// @cite astrom2006 -- Ch. 3.3 (integral action, discretisation methods)
    auto compute_integral_term(const vector_t& e, Scalar dt) -> std::pair<vector_t, vector_t>
    {
        vector_t integral_increment = vector_t::Zero();
        if(!m_integral_frozen)
        {
            if constexpr(detail::contains_v<forward_euler, Policies...>)
                integral_increment = (m_ki.cwiseProduct(m_prev_error) * dt).eval();
            else if constexpr(detail::contains_v<tustin, Policies...>)
            {
                auto avg = ((e + m_prev_error) * Scalar{0.5}).eval();
                integral_increment = (m_ki.cwiseProduct(avg) * dt).eval();
            }
            else
                integral_increment = (m_ki.cwiseProduct(e) * dt).eval();
            m_integral = (m_integral + integral_increment).eval();
        }
        return {integral_increment, m_integral};
    }

    /// @cite astrom2006 -- Ch. 3.4 (derivative action, setpoint weighting c parameter)
    auto compute_derivative_term(const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        vector_t d = vector_t::Zero();
        if(!m_first_step)
        {
            d = compute_raw_derivative(filtered_sp, filtered_meas, dt);
            d = apply_derivative_filter(d, dt);
        }
        else if constexpr(detail::contains_v<deriv_filter, Policies...>)
            m_prev_deriv_filtered = d;
        return d;
    }

    auto compute_raw_derivative(const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        if(m_cfg.derivative_on_error)
        {
            vector_t ed_curr, ed_prev;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                ed_curr[i] = m_cfg.c[i] * filtered_sp[i] - filtered_meas[i];
                ed_prev[i] = m_cfg.c[i] * m_prev_sp[i] - m_prev_meas[i];
            }
            return m_kd.cwiseProduct((ed_curr - ed_prev) / dt).eval();
        }
        else
        {
            auto dm_dt = ((filtered_meas - m_prev_meas) / dt).eval();
            return (-m_kd.cwiseProduct(dm_dt)).eval();
        }
    }

    auto apply_derivative_filter(vector_t d, Scalar dt) -> vector_t
    {
        if constexpr(detail::contains_v<deriv_filter, Policies...>)
        {
            const auto& df_cfg = m_cfg.template policy<deriv_filter>();
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                if(df_cfg.n[static_cast<std::size_t>(i)] > Scalar{0} && m_kp[i] != Scalar{0})
                {
                    Scalar tf = m_kd[i] / (m_kp[i] * df_cfg.n[static_cast<std::size_t>(i)]);
                    Scalar alpha = tf / (tf + dt);
                    d[i] = alpha * m_prev_deriv_filtered[i] + (Scalar{1} - alpha) * d[i];
                }
            m_prev_deriv_filtered = d;
        }
        return d;
    }

    auto compute_raw_output(const vector_t& p, const vector_t& integral, const vector_t& d, const vector_t& sp, Scalar dt) -> vector_t
    {
        return apply_feed_forward((p + integral + d).eval(), sp, dt);
    }

    auto apply_feed_forward(vector_t u_raw, const vector_t& sp, Scalar dt) -> vector_t
    {
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                auto ff = m_cfg.template policy<ff_policy_t>().ff_func(sp, dt);
                return (u_raw + ff).eval();
            }
        }
        return u_raw;
    }

    auto apply_rate_limit(vector_t u_raw, Scalar dt) -> vector_t
    {
        if constexpr(detail::contains_v<rate_limit, Policies...>)
        {
            auto delta = (u_raw - m_prev_output).eval();
            const auto& rl_cfg = m_cfg.template policy<rate_limit>();
            vector_t lo, hi;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(rl_cfg.rate_max[static_cast<std::size_t>(i)] > Scalar{0} && rl_cfg.rate_max[static_cast<std::size_t>(i)] != std::numeric_limits<Scalar>::infinity())
                {
                    hi[i] = rl_cfg.rate_max[static_cast<std::size_t>(i)] * dt;
                    lo[i] = -hi[i];
                }
                else
                {
                    hi[i] = std::numeric_limits<Scalar>::max();
                    lo[i] = std::numeric_limits<Scalar>::lowest();
                }
            }
            return (m_prev_output + delta.cwiseMax(lo).cwiseMin(hi)).eval();
        }
        return u_raw;
    }

    void apply_anti_windup(const vector_t& u_sat, const vector_t& u_unconstrained, const vector_t& e, const vector_t& integral_increment, Scalar dt)
    {
        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
        {
            using AW = detail::find_policy_t<anti_windup, Policies...>;
            if constexpr(std::is_same_v<AW, anti_windup<back_calc>>)
            {
                auto feedback = kb_.cwiseProduct((u_sat - u_unconstrained).eval()).eval();
                m_integral = (m_integral + feedback * dt).eval();
            }
            else if constexpr(std::is_same_v<AW, anti_windup<clamping>>)
            {
                // Gate the integral undo per channel on that channel's own constraint
                // (applied output differs from the unconstrained command), so a
                // saturated channel does not freeze an unsaturated one.
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                    if(u_sat[i] != u_unconstrained[i]
                        && ((e[i] > Scalar{0} && m_integral[i] > Scalar{0}) || (e[i] < Scalar{0} && m_integral[i] < Scalar{0})))
                        m_integral[i] -= integral_increment[i];
            }
            else if constexpr(std::is_same_v<AW, anti_windup<conditional_integration>>)
            {
                const auto& ci_cfg = m_cfg.template policy<AW>();
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    if(abs_e > ci_cfg.error_threshold[static_cast<std::size_t>(i)])
                        m_integral[i] -= integral_increment[i];
                }
            }
        }
    }

    void update_state(const vector_t& e, const vector_t& filtered_meas, const vector_t& filtered_sp, const vector_t& output)
    {
        m_prev_prev_error = m_prev_error;
        m_prev_error = e;
        m_prev_meas = filtered_meas;
        m_prev_sp = filtered_sp;
        m_prev_output = output;
        m_first_step = false;
        m_perf.set_first_step(false);
    }

    config_type m_cfg;
    vector_t m_kp = vector_t::Zero();
    vector_t m_ki = vector_t::Zero();
    vector_t m_kd = vector_t::Zero();
    vector_t kb_ = vector_t::Zero();
    vector_t m_integral = vector_t::Zero();
    vector_t m_prev_error = vector_t::Zero();
    vector_t m_prev_prev_error = vector_t::Zero();
    vector_t m_prev_meas = vector_t::Zero();
    vector_t m_prev_sp = vector_t::Zero();
    vector_t m_prev_output = vector_t::Zero();
    vector_t m_accumulated_output = vector_t::Zero();
    vector_t m_prev_ff = vector_t::Zero();
    vector_t m_filtered_sp = vector_t::Zero();
    vector_t m_filtered_meas = vector_t::Zero();
    vector_t m_prev_deriv_filtered = vector_t::Zero();
    pid_performance_tracker<Scalar, NY, Policies...> m_perf;
    bool m_first_step{true};
    bool m_integral_frozen{false};
    bool m_saturated{false};
};

}

#endif
