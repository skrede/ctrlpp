#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_CYCLE_STATE_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_CYCLE_STATE_H

#include "ctrlpp/control/detail/pid_gains.h"
#include "ctrlpp/control/detail/pid_integrator.h"
#include "ctrlpp/control/detail/pid_output_stage.h"
#include "ctrlpp/control/detail/pid_input_filters.h"
#include "ctrlpp/control/detail/pid_derivative_path.h"

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"
#include "ctrlpp/control/pid_performance.h"

#include "ctrlpp/types.h"

#include <cstddef>
#include <utility>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_cycle_state
{
    using vector_t = Vector<Scalar, NY>;
    using gains_type = pid_gains<Scalar, NY>;
    using config_type = pid_config<Scalar, NY, Policies...>;
    using tracker_type = pid_performance_tracker<Scalar, NY, Policies...>;

public:
    pid_cycle_state()
        : m_prev_sp{vector_t::Zero()}
        , m_prev_meas{vector_t::Zero()}
        , m_prev_error{vector_t::Zero()}
        , m_first_step{true}
        , m_saturated{false}
    {
    }

    const vector_t& error() const { return m_prev_error; }
    const vector_t& integral() const { return m_integrator.integral(); }
    bool saturated() const { return m_saturated; }
    tracker_type& performance() { return m_performance; }
    const tracker_type& performance() const { return m_performance; }

    void assign_integral(const vector_t& value) { m_integrator.assign(value); }

    vector_t advance(config_type& cfg, const gains_type& gains, bool frozen, const vector_t& sp, const vector_t& meas, Scalar dt, vector_t& last_non_integral)
    {
        auto filtered_sp = m_filters.apply_setpoint_filter(cfg, sp, dt);
        auto filtered_meas = m_filters.apply_pv_filter(cfg, meas, dt);
        auto e = (filtered_sp - filtered_meas).eval();
        m_performance.accumulate(e, dt, cfg);
        if constexpr(contains_v<velocity_form, Policies...>)
            return compute_velocity_form(cfg, gains, frozen, e, sp, filtered_sp, filtered_meas, dt);
        else
            return compute_position_form(cfg, gains, frozen, e, sp, filtered_sp, filtered_meas, dt, last_non_integral);
    }

    bool track(const vector_t& target, const vector_t& last_non_integral)
    {
        if constexpr(!contains_v<velocity_form, Policies...>)
        {
            auto next_integral = (target - last_non_integral).eval();
            if(!next_integral.allFinite())
                return false;
            m_integrator.track(std::move(next_integral));
        }
        m_output.track(target);
        return true;
    }

    // Each component gates its members on the policy that makes them live, so a value seeded into a member the
    // composed controller never reads cannot cause a rejection the output would not have suffered from.
    bool finite() const
    {
        if(!m_prev_error.allFinite() || !m_prev_meas.allFinite() || !m_prev_sp.allFinite())
            return false;
        return m_integrator.finite() && m_output.finite() && m_derivative.finite() && m_filters.finite();
    }

private:
    vector_t m_prev_sp;
    vector_t m_prev_meas;
    vector_t m_prev_error;
    tracker_type m_performance;
    pid_output_stage<Scalar, NY, Policies...> m_output;
    pid_integrator<Scalar, NY, Policies...> m_integrator;
    pid_input_filters<Scalar, NY, Policies...> m_filters;
    pid_derivative_path<Scalar, NY, Policies...> m_derivative;
    // The flags come last: declared first, they cost the float cycle a register save on Cortex-M7.
    bool m_first_step;
    bool m_saturated;

    // The incremental form as the first difference of the positional terms, so step 1 has no
    // derivative increment and each cycle's own dt divides its own difference.
    vector_t compute_velocity_form(
            config_type& cfg, const gains_type& gains, bool frozen, const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt)
    {
        auto dp = (compute_proportional_term(cfg, gains, filtered_sp, filtered_meas) - compute_proportional_term(cfg, gains, m_prev_sp, m_prev_meas)).eval();
        auto d = m_derivative.compute_derivative_term(cfg, gains, filtered_sp, filtered_meas, m_prev_sp, m_prev_meas, m_first_step, dt);
        auto delta_u = (dp + m_integrator.integral_increment(gains, frozen, m_prev_error, e, dt) + m_derivative.increment(d)).eval();
        auto with_ff = m_output.apply_feed_forward_velocity(cfg, delta_u, sp, dt);
        auto emitted = m_output.accumulate_output(cfg, m_output.clip_increment(cfg, with_ff, dt).value, m_saturated);
        update_state(e, filtered_meas, filtered_sp, emitted);
        return emitted;
    }

    // @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3-4
    vector_t compute_position_form(config_type& cfg, const gains_type& gains, bool frozen, const vector_t& e, const vector_t& sp, const vector_t& filtered_sp,
            const vector_t& filtered_meas, Scalar dt, vector_t& last_non_integral)
    {
        auto p = compute_proportional_term(cfg, gains, filtered_sp, filtered_meas);
        auto integral = m_integrator.compute_integral_term(cfg, gains, frozen, m_prev_error, e, dt);
        auto d = m_derivative.compute_derivative_term(cfg, gains, filtered_sp, filtered_meas, m_prev_sp, m_prev_meas, m_first_step, dt);
        auto u_unconstrained = compute_raw_output(cfg, p, integral, d, sp, dt, last_non_integral);
        auto limited = m_output.apply_rate_limit(cfg, u_unconstrained, dt);
        auto u_sat = limited.value.cwiseMax(cfg.output_min).cwiseMin(cfg.output_max).eval();
        auto saturated = (limited.value.array() < cfg.output_min.array() || limited.value.array() > cfg.output_max.array()).eval();
        m_saturated = saturated.any();
        // Anti-windup feeds back against the command taken before the rate limiter, so a rate-limited ramp,
        // where the applied output already equals the limited command, engages it as well as saturation does.
        m_integrator.apply_anti_windup(gains, frozen, u_sat, u_unconstrained, (limited.clipped || saturated).eval(), dt);
        update_state(e, filtered_meas, filtered_sp, u_sat);
        return u_sat;
    }

    // @cite astrom2006 -- Ch. 3.5 (setpoint weighting b parameter)
    static vector_t compute_proportional_term(const config_type& cfg, const gains_type& gains, const vector_t& filtered_sp, const vector_t& filtered_meas)
    {
        vector_t ep;
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            ep[i] = cfg.b[i] * filtered_sp[i] - filtered_meas[i];
        return gains.kp.cwiseProduct(ep).eval();
    }

    static vector_t compute_raw_output(config_type& cfg, const vector_t& p, const vector_t& integral, const vector_t& d, const vector_t& sp, Scalar dt, vector_t& last_non_integral)
    {
        auto u_raw = (p + integral + d).eval();
        last_non_integral = (p + d).eval();
        if constexpr(has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                // The feed-forward callable may be stateful, so it is evaluated once per cycle.
                const vector_t ff = cfg.template policy<ff_policy_t>().feed_forward_cb(sp, dt);
                last_non_integral += ff;
                return (u_raw + ff).eval();
            }
        }
        return u_raw;
    }

    void update_state(const vector_t& e, const vector_t& filtered_meas, const vector_t& filtered_sp, const vector_t& output)
    {
        m_prev_error = e;
        m_prev_meas = filtered_meas;
        m_prev_sp = filtered_sp;
        m_output.record(output);
        m_first_step = false;
    }
};

}

#endif
