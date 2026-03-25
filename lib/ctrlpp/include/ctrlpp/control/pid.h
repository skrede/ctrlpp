#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

/// @brief Policy-based PID controller with compile-time feature composition.
///
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/types.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename... Policies>
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
    }

    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt) -> vector_t
    {
        if(dt <= Scalar{0})
            return m_prev_output;

        auto filtered_sp = apply_setpoint_filter(sp, dt);
        auto filtered_meas = apply_pv_filter(meas, dt);
        auto e = (filtered_sp - filtered_meas).eval();

        accumulate_performance_metrics(e, dt);

        if constexpr(detail::contains_v<velocity_form, Policies...>)
            return compute_velocity_form(e, sp, filtered_sp, filtered_meas, dt);
        else
            return compute_position_form(e, sp, filtered_sp, filtered_meas, dt);
    }

    vector_t compute(const vector_t& sp, const vector_t& meas, Scalar dt, const vector_t& tracking_signal)
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
        vector_t ki_old = m_ki;
        m_cfg = new_cfg;
        compute_internal_gains(new_cfg);
        rescale_integral_bumpless(ki_old);

        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
            initialize_back_calc_gains();
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
            reset_metrics();
    }

    void freeze_integral(bool freeze = true) { m_integral_frozen = freeze; }

    void set_integral(const vector_t& val) { m_integral = val; }

    template <typename Metric>
    const vector_t& metric() const
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

    bool oscillating() const
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

    void reset_metrics()
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

private:
    // --- Gain initialization sub-steps ---

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
                    if(m_kd[i] != Scalar{0})
                        kb_[i] = std::sqrt(m_ki[i] * m_kd[i]);
                    else
                        kb_[i] = m_ki[i];
                }
            }
        }
    }

    void rescale_integral_bumpless(const vector_t& ki_old)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if(m_ki[i] != Scalar{0} && ki_old[i] != Scalar{0})
                m_integral[i] = m_integral[i] * ki_old[i] / m_ki[i];
            else if(m_ki[i] == Scalar{0})
                m_integral[i] = Scalar{0};
        }
    }

    // --- Input filtering sub-steps ---

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

    // --- Performance metric accumulation ---

    void accumulate_performance_metrics(const vector_t& e, Scalar dt)
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

    // --- Velocity form computation ---

    auto compute_velocity_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        auto dp = m_kp.cwiseProduct(e - m_prev_error).eval();
        auto di = (m_ki.cwiseProduct(e) * dt).eval();
        auto d_num = (e - m_prev_error * Scalar{2} + m_prev_prev_error).eval();
        auto dd = m_kd.cwiseProduct(d_num / dt).eval();
        auto delta_u = (dp + di + dd).eval();

        delta_u = apply_feed_forward_velocity(delta_u, sp, dt);
        delta_u = delta_u.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();

        update_state(e, filtered_meas, filtered_sp, delta_u);
        return delta_u;
    }

    auto apply_feed_forward_velocity(vector_t delta_u, const vector_t& sp, Scalar dt) -> vector_t
    {
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                auto ff = m_cfg.template policy<ff_policy_t>().ff_func(sp, dt);
                return (delta_u + ff).eval();
            }
        }
        return delta_u;
    }

    // --- Position form computation ---
    /// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3-4

    auto compute_position_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        auto p = compute_proportional_term(filtered_sp, filtered_meas);
        auto [integral_increment, updated_integral] = compute_integral_term(e, dt);
        auto d = compute_derivative_term(filtered_sp, filtered_meas, dt);
        auto u_raw = compute_raw_output(p, updated_integral, d, sp, dt);
        u_raw = apply_rate_limit(u_raw, dt);

        auto u_sat = u_raw.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        m_saturated = (u_sat.array() != u_raw.array()).any();

        apply_anti_windup(u_sat, u_raw, e, integral_increment, dt);
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
            vector_t ed_curr;
            vector_t ed_prev;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                ed_curr[i] = m_cfg.c[i] * filtered_sp[i] - filtered_meas[i];
                ed_prev[i] = m_cfg.c[i] * m_prev_sp[i] - m_prev_meas[i];
            }
            auto de = (ed_curr - ed_prev).eval();
            return m_kd.cwiseProduct(de / dt).eval();
        }
        else
        {
            auto dm = (filtered_meas - m_prev_meas).eval();
            auto dm_dt = (dm / dt).eval();
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
        auto u_raw = (p + integral + d).eval();
        return apply_feed_forward(u_raw, sp, dt);
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
            auto [lo, hi] = compute_rate_delta_bounds(dt);
            auto clamped_delta = delta.cwiseMax(lo).cwiseMin(hi).eval();
            return (m_prev_output + clamped_delta).eval();
        }
        return u_raw;
    }

    auto compute_rate_delta_bounds(Scalar dt) const -> std::pair<vector_t, vector_t>
    {
        const auto& rl_cfg = m_cfg.template policy<rate_limit>();
        vector_t max_delta;
        vector_t neg_max_delta;
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if(rl_cfg.rate_max[static_cast<std::size_t>(i)] > Scalar{0} && rl_cfg.rate_max[static_cast<std::size_t>(i)] != std::numeric_limits<Scalar>::infinity())
            {
                max_delta[i] = rl_cfg.rate_max[static_cast<std::size_t>(i)] * dt;
                neg_max_delta[i] = -max_delta[i];
            }
            else
            {
                max_delta[i] = std::numeric_limits<Scalar>::max();
                neg_max_delta[i] = std::numeric_limits<Scalar>::lowest();
            }
        }
        return {neg_max_delta, max_delta};
    }

    void apply_anti_windup(const vector_t& u_sat, const vector_t& u_raw, const vector_t& e, const vector_t& integral_increment, Scalar dt)
    {
        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
        {
            using AW = detail::find_policy_t<anti_windup, Policies...>;
            if constexpr(std::is_same_v<AW, anti_windup<back_calc>>)
                apply_back_calc_windup(u_sat, u_raw, dt);
            else if constexpr(std::is_same_v<AW, anti_windup<clamping>>)
                apply_clamping_windup(e, integral_increment);
            else if constexpr(std::is_same_v<AW, anti_windup<conditional_integration>>)
                apply_conditional_integration_windup(e, integral_increment);
        }
    }

    /// @cite astrom2006 -- Ch. 3.5 (back-calculation anti-windup)
    void apply_back_calc_windup(const vector_t& u_sat, const vector_t& u_raw, Scalar dt)
    {
        auto sat_diff = (u_sat - u_raw).eval();
        auto feedback = kb_.cwiseProduct(sat_diff).eval();
        m_integral = (m_integral + feedback * dt).eval();
    }

    /// @cite astrom2006 -- Ch. 3.5 (clamping anti-windup)
    void apply_clamping_windup(const vector_t& e, const vector_t& integral_increment)
    {
        if(m_saturated)
        {
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if((e[i] > Scalar{0} && m_integral[i] > Scalar{0}) || (e[i] < Scalar{0} && m_integral[i] < Scalar{0}))
                    m_integral[i] -= integral_increment[i];
            }
        }
    }

    /// @cite astrom2006 -- Ch. 3.5 (conditional integration anti-windup)
    void apply_conditional_integration_windup(const vector_t& e, const vector_t& integral_increment)
    {
        using AW = detail::find_policy_t<anti_windup, Policies...>;
        const auto& ci_cfg = m_cfg.template policy<AW>();
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
            if(abs_e > ci_cfg.error_threshold[static_cast<std::size_t>(i)])
                m_integral[i] -= integral_increment[i];
        }
    }

    // --- State update ---

    void update_state(const vector_t& e, const vector_t& filtered_meas, const vector_t& filtered_sp, const vector_t& output)
    {
        m_prev_prev_error = m_prev_error;
        m_prev_error = e;
        m_prev_meas = filtered_meas;
        m_prev_sp = filtered_sp;
        m_prev_output = output;
        m_first_step = false;
    }

    // --- Member data ---

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
    vector_t m_filtered_sp = vector_t::Zero();
    vector_t m_filtered_meas = vector_t::Zero();
    vector_t m_prev_deriv_filtered = vector_t::Zero();
    vector_t m_iae = vector_t::Zero();
    vector_t m_ise = vector_t::Zero();
    vector_t m_itae = vector_t::Zero();
    vector_t m_zero_crossings = vector_t::Zero();
    vector_t m_prev_error_sign = vector_t::Zero();
    Scalar m_accumulated_time{0};
    Scalar m_osc_threshold{5.0};
    bool m_first_step{true};
    bool m_integral_frozen{false};
    bool m_saturated{false};
};

}

#endif
