#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

#include "ctrlpp/types.h"
#include "ctrlpp/pid_config.h"
#include "ctrlpp/pid_policies.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp {

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename... Policies>
class pid
{
public:
    using config_type = pid_config<Scalar, NY, Policies...>;
    using vector_t = Vector<Scalar, NY>;

    explicit pid(const config_type &cfg)
        : m_cfg{cfg}
    {
        if constexpr(detail::contains_v<isa_form, Policies...>)
        {
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(m_cfg.ki[i] != Scalar{0})
                    m_ki[i] = m_cfg.kp[i] / m_cfg.ki[i];
                else
                    m_ki[i] = Scalar{0};
                m_kd[i] = m_cfg.kp[i] * m_cfg.kd[i];
                m_kp[i] = m_cfg.kp[i];
            }
        }
        else
        {
            m_kp = m_cfg.kp;
            m_ki = m_cfg.ki;
            m_kd = m_cfg.kd;
        }

        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
        {
            using AW = detail::find_policy_t<anti_windup, Policies...>;
            if constexpr(std::is_same_v<AW, anti_windup<back_calc>>)
            {
                const auto &aw_cfg = m_cfg.template policy<AW>();
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
    }

    auto compute(const vector_t &sp, const vector_t &meas, Scalar dt) -> vector_t
    {
        if(dt <= Scalar{0})
            return m_prev_output;

        // Input filtering: setpoint filter
        auto filtered_sp = sp;
        if constexpr(detail::contains_v<setpoint_filter, Policies...>)
        {
            const auto &tf_sp = m_cfg.template policy<setpoint_filter>().tf;
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
            filtered_sp = m_filtered_sp;
        }

        // Input filtering: process variable filter
        auto filtered_meas = meas;
        if constexpr(detail::contains_v<pv_filter, Policies...>)
        {
            const auto &tf_pv = m_cfg.template policy<pv_filter>().tf;
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
            filtered_meas = m_filtered_meas;
        }

        // Error (full, using filtered signals)
        auto e = (filtered_sp - filtered_meas).eval();

        // Performance metric accumulation
        if constexpr(detail::has_policy_v<perf_assessment, Policies...>)
        {
            using PA = detail::find_policy_t<perf_assessment, Policies...>;
            m_accumulated_time += dt;

            if constexpr(detail::perf_has_metric_v<IAE, PA>)
            {
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    m_iae[i] += abs_e * dt;
                }
            }
            if constexpr(detail::perf_has_metric_v<ISE, PA>)
            {
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                    m_ise[i] += e[i] * e[i] * dt;
            }
            if constexpr(detail::perf_has_metric_v<ITAE, PA>)
            {
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    m_itae[i] += m_accumulated_time * abs_e * dt;
                }
            }
            if constexpr(detail::perf_has_metric_v<oscillation_detect, PA>)
            {
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                {
                    Scalar sign_e = (e[i] > Scalar{0})
                                        ? Scalar{1}
                                        : (e[i] < Scalar{0})
                                        ? Scalar{-1}
                                        : Scalar{0};
                    if(!m_first_step && sign_e != Scalar{0} && m_prev_error_sign[i] != Scalar{0} && sign_e != m_prev_error_sign[i])
                        m_zero_crossings[i] += Scalar{1};
                    if(sign_e != Scalar{0})
                        m_prev_error_sign[i] = sign_e;
                }
            }
        }

        if constexpr(detail::contains_v<velocity_form, Policies...>)
        {
            // Velocity (incremental) form: outputs delta_u per step
            // P contribution: Kp * (e(k) - e(k-1))
            auto dp = m_kp.cwiseProduct(e - m_prev_error).eval();

            // I contribution: Ki * e(k) * dt (no accumulator — naturally anti-windup)
            auto di = (m_ki.cwiseProduct(e) * dt).eval();

            // D contribution: Kd * (e(k) - 2*e(k-1) + e(k-2)) / dt
            auto d_num = (e - m_prev_error * Scalar{2} + m_prev_prev_error).eval();
            auto dd = m_kd.cwiseProduct(d_num / dt).eval();

            auto delta_u = (dp + di + dd).eval();

            // Feed-forward (adds to delta_u)
            if constexpr(detail::has_policy_v<feed_forward, Policies...>)
            {
                using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
                if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
                {
                    auto ff = m_cfg.template policy<ff_policy_t>().ff_func(sp, dt);
                    delta_u = (delta_u + ff).eval();
                }
            }

            // Output clamp on delta_u
            delta_u = delta_u.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();

            // Update state
            m_prev_prev_error = m_prev_error;
            m_prev_error = e;
            m_prev_meas = filtered_meas;
            m_prev_sp = filtered_sp;
            m_prev_output = delta_u;
            m_first_step = false;

            return delta_u;
        }
        else
        {
            // Position (absolute) form
            // P term with setpoint weighting: ep = b*sp - meas
            vector_t ep;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                ep[i] = m_cfg.b[i] * filtered_sp[i] - filtered_meas[i];
            auto p = m_kp.cwiseProduct(ep).eval();

            // I term -- always uses full error (sp - meas)
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

            // D term with optional setpoint weighting (c parameter) and derivative filter
            vector_t d = vector_t::Zero();
            if(!m_first_step)
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
                    d = m_kd.cwiseProduct(de / dt).eval();
                }
                else
                {
                    auto dm = (filtered_meas - m_prev_meas).eval();
                    auto dm_dt = (dm / dt).eval();
                    d = (-m_kd.cwiseProduct(dm_dt)).eval();
                }

                // Derivative filter (first-order low-pass)
                if constexpr(detail::contains_v<deriv_filter, Policies...>)
                {
                    const auto &df_cfg = m_cfg.template policy<deriv_filter>();
                    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                        if(df_cfg.n[static_cast<std::size_t>(i)] > Scalar{0} && m_kp[i] != Scalar{0})
                        {
                            Scalar tf = m_kd[i] / (m_kp[i] * df_cfg.n[static_cast<std::size_t>(i)]);
                            Scalar alpha = tf / (tf + dt);
                            d[i] = alpha * m_prev_deriv_filtered[i] + (Scalar{1} - alpha) * d[i];
                        }
                    m_prev_deriv_filtered = d;
                }
            }
            else if constexpr(detail::contains_v<deriv_filter, Policies...>)
                m_prev_deriv_filtered = d;

            // Sum
            auto u_raw = (p + m_integral + d).eval();

            // Feed-forward
            if constexpr(detail::has_policy_v<feed_forward, Policies...>)
            {
                using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
                if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
                {
                    auto ff = m_cfg.template policy<ff_policy_t>().ff_func(sp, dt);
                    u_raw = (u_raw + ff).eval();
                }
            }

            // Rate limit (if rate_limit policy present)
            if constexpr(detail::contains_v<rate_limit, Policies...>)
            {
                const auto &rl_cfg = m_cfg.template policy<rate_limit>();
                auto delta = (u_raw - m_prev_output).eval();
                vector_t max_delta;
                vector_t neg_max_delta;
                for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                {
                    if(rl_cfg.rate_max[static_cast<std::size_t>(i)] > Scalar{0} &&
                        rl_cfg.rate_max[static_cast<std::size_t>(i)] != std::numeric_limits<Scalar>::infinity())
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
                auto clamped_delta = delta.cwiseMax(neg_max_delta).cwiseMin(max_delta).eval();
                u_raw = (m_prev_output + clamped_delta).eval();
            }

            // Output clamp
            auto u_sat = u_raw.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();

            // Saturated flag
            m_saturated = (u_sat.array() != u_raw.array()).any();

            // Anti-windup feedback (adjusts integral for next step)
            if constexpr(detail::has_policy_v<anti_windup, Policies...>)
            {
                using AW = detail::find_policy_t<anti_windup, Policies...>;
                if constexpr(std::is_same_v<AW, anti_windup<back_calc>>)
                {
                    auto sat_diff = (u_sat - u_raw).eval();
                    auto feedback = kb_.cwiseProduct(sat_diff).eval();
                    m_integral = (m_integral + feedback * dt).eval();
                }
                else if constexpr(std::is_same_v<AW, anti_windup<clamping>>)
                {
                    if(m_saturated)
                    {
                        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                        {
                            if((e[i] > Scalar{0} && m_integral[i] > Scalar{0}) ||
                                (e[i] < Scalar{0} && m_integral[i] < Scalar{0}))
                            {
                                m_integral[i] -= integral_increment[i];
                            }
                        }
                    }
                }
                else if constexpr(std::is_same_v<AW, anti_windup<conditional_integration>>)
                {
                    const auto &ci_cfg = m_cfg.template policy<AW>();
                    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                    {
                        Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                        if(abs_e > ci_cfg.error_threshold[static_cast<std::size_t>(i)])
                        {
                            m_integral[i] -= integral_increment[i];
                        }
                    }
                }
            }

            // Update state
            m_prev_error = e;
            m_prev_meas = filtered_meas;
            m_prev_sp = filtered_sp;
            m_prev_output = u_sat;
            m_first_step = false;

            return u_sat;
        }
    }

    vector_t compute(const vector_t &sp, const vector_t &meas, Scalar dt, const vector_t &tracking_signal)
    {
        // Compute normal output first
        auto u = compute(sp, meas, dt);

        if(dt <= Scalar{0})
            return u;

        if constexpr(!detail::contains_v<velocity_form, Policies...>)
        {
            // Tracking: set integral so that output matches tracking_signal
            auto non_integral = (u - m_integral).eval();
            m_integral = (tracking_signal - non_integral).eval();
        }

        return u;
    }

    void set_params(const config_type &new_cfg)
    {
        // Bumpless integral rescaling: integral_new = integral_old * ki_old / ki_new
        vector_t ki_old = m_ki;

        // Update config
        m_cfg = new_cfg;

        // Recompute internal gains
        if constexpr(detail::contains_v<isa_form, Policies...>)
        {
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                if(new_cfg.ki[i] != Scalar{0})
                    m_ki[i] = new_cfg.kp[i] / new_cfg.ki[i];
                else
                    m_ki[i] = Scalar{0};
                m_kd[i] = new_cfg.kp[i] * new_cfg.kd[i];
                m_kp[i] = new_cfg.kp[i];
            }
        }
        else
        {
            m_kp = new_cfg.kp;
            m_ki = new_cfg.ki;
            m_kd = new_cfg.kd;
        }

        // Rescale integral for bumpless gain change
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if(m_ki[i] != Scalar{0} && ki_old[i] != Scalar{0})
                m_integral[i] = m_integral[i] * ki_old[i] / m_ki[i];
            else if(m_ki[i] == Scalar{0})
                m_integral[i] = Scalar{0};
        }

        // Recompute back_calc Kb if applicable
        if constexpr(detail::has_policy_v<anti_windup, Policies...>)
        {
            using AW = detail::find_policy_t<anti_windup, Policies...>;
            if constexpr(std::is_same_v<AW, anti_windup<back_calc>>)
            {
                const auto &aw_cfg = m_cfg.template policy<AW>();
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
    }

    const vector_t &error() const
    {
        return m_prev_error;
    }

    const vector_t &integral() const
    {
        return m_integral;
    }

    const config_type &params() const
    {
        return m_cfg;
    }

    bool saturated() const
    {
        return m_saturated;
    }

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

    void freeze_integral(bool freeze = true)
    {
        m_integral_frozen = freeze;
    }

    void set_integral(const vector_t &val)
    {
        m_integral = val;
    }

    template <typename Metric>
    const vector_t &metric() const requires detail::has_policy_v<perf_assessment, Policies...>
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

    bool oscillating() const requires detail::has_policy_v<perf_assessment, Policies...>
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

    void reset_metrics() requires detail::has_policy_v<perf_assessment, Policies...>
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
