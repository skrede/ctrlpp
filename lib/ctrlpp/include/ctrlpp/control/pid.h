#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/control/pid_error.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"
#include "ctrlpp/control/pid_performance.h"

#include "ctrlpp/control/detail/pid_gains.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <cstdint>

namespace ctrlpp
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid
{
    static_assert(!(detail::contains_v<velocity_form, Policies...> && detail::has_policy_v<anti_windup, Policies...>),
            "anti_windup does not compose with velocity_form: clamping the accumulated output already prevents windup");

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    using mask_t = Eigen::Array<bool, static_cast<int>(NY), 1>;
    using direction_t = Eigen::Array<std::int8_t, static_cast<int>(NY), 1>;

    struct limited_command
    {
        Vector<Scalar, NY> value;
        mask_t clipped;
    };

    struct cycle_state
    {
        Vector<Scalar, NY> integral;
        Vector<Scalar, NY> prev_error;
        Vector<Scalar, NY> prev_derivative;
        Vector<Scalar, NY> prev_meas;
        Vector<Scalar, NY> prev_sp;
        Vector<Scalar, NY> prev_output;
        Vector<Scalar, NY> accumulated_output;
        Vector<Scalar, NY> prev_ff;
        Vector<Scalar, NY> filtered_sp;
        Vector<Scalar, NY> filtered_meas;
        Vector<Scalar, NY> prev_deriv_filtered;
        direction_t constraint_direction;
        pid_performance_tracker<Scalar, NY, Policies...> performance;
        bool first_step;
        bool saturated;
    };

public:
    using config_type = pid_config<Scalar, NY, Policies...>;
    using vector_t = Vector<Scalar, NY>;

    static expected<pid, pid_config_error> create(const config_type& cfg)
    {
        auto derived = detail::derive_pid_gains(cfg);
        if(!derived)
            return unexpected(derived.error());
        return pid{unchecked_t{}, cfg, *derived};
    }

    expected<vector_t, pid_step_error> compute(const vector_t& sp, const vector_t& meas, Scalar dt)
    {
        if(const auto step = check_step(sp, meas, dt); !step)
            return unexpected(latch_health(step.error()));

        auto const previous = capture_cycle_state();
        auto filtered_sp = apply_setpoint_filter(sp, dt);
        auto filtered_meas = apply_pv_filter(meas, dt);
        auto e = (filtered_sp - filtered_meas).eval();
        m_perf.accumulate(e, dt, m_cfg);
        auto result = compute_form(e, sp, filtered_sp, filtered_meas, dt);
        if(!result.allFinite() || !carried_state_finite() || !m_perf.all_finite())
        {
            restore_cycle_state(previous);
            return unexpected(pid_step_error::non_finite_result);
        }
        return result;
    }

    // The tracking signal is checked before the cycle runs because a successful cycle writes it straight into
    // carried state, after the cycle has committed the rest.
    expected<vector_t, pid_step_error> compute(const vector_t& sp, const vector_t& meas, Scalar dt, const vector_t& tracking_signal)
    {
        if(!tracking_signal.allFinite())
            return unexpected(pid_step_error::non_finite_tracking_signal);

        auto const previous = capture_cycle_state();
        auto u = compute(sp, meas, dt);
        if(!u)
            return u;
        if(!track(tracking_signal.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval()))
        {
            restore_cycle_state(previous);
            return unexpected(pid_step_error::non_finite_result);
        }
        return u;
    }

    expected<void, pid_config_error> set_params(const config_type& new_cfg)
    {
        auto derived = detail::derive_pid_gains(new_cfg);
        if(!derived)
            return unexpected(derived.error());
        // The integral is stored in output units, so its contribution is already continuous across a
        // gain change: leaving it untouched is the bumpless transfer, rescaling it would inject a bump.
        m_cfg = new_cfg;
        m_gains = *derived;
        return {};
    }

    const vector_t& error() const { return m_prev_error; }
    const vector_t& integral() const { return m_integral; }
    const config_type& params() const { return m_cfg; }
    bool saturated() const { return m_saturated; }

    pid_health health() const { return m_health; }

    void reset()
    {
        m_integral = vector_t::Zero();
        m_prev_error = vector_t::Zero();
        m_prev_meas = vector_t::Zero();
        m_prev_sp = vector_t::Zero();
        m_prev_output = vector_t::Zero();
        m_accumulated_output = vector_t::Zero();
        m_prev_ff = vector_t::Zero();
        m_prev_derivative = vector_t::Zero();
        m_filtered_sp = vector_t::Zero();
        m_filtered_meas = vector_t::Zero();
        m_prev_deriv_filtered = vector_t::Zero();
        m_first_step = true;
        m_integral_frozen = false;
        m_saturated = false;
        m_constraint_direction = direction_t::Zero();
        m_health = pid_health::ok;
        if constexpr(detail::has_policy_v<perf_assessment, Policies...>)
            m_perf.reset();
    }

    void freeze_integral(bool freeze = true) { m_integral_frozen = freeze; }

    expected<void, pid_config_error> set_integral(const vector_t& val)
        requires(!detail::contains_v<velocity_form, Policies...>)
    {
        if(!val.allFinite())
            return unexpected(pid_config_error::non_finite_integral);
        m_integral = val;
        return {};
    }

    template <typename Metric>
    const vector_t& metric() const
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_perf.template metric<Metric>();
    }

    bool oscillating() const
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_perf.oscillating(m_cfg);
    }

    void reset_metrics()
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        m_perf.reset();
    }

private:
    config_type m_cfg;
    detail::pid_gains<Scalar, NY> m_gains;
    vector_t m_integral;
    vector_t m_prev_error;
    vector_t m_prev_derivative;
    vector_t m_prev_meas;
    vector_t m_prev_sp;
    vector_t m_prev_output;
    vector_t m_accumulated_output;
    vector_t m_prev_ff;
    vector_t m_filtered_sp;
    vector_t m_filtered_meas;
    vector_t m_prev_deriv_filtered;
    direction_t m_constraint_direction;
    // Per-cycle scratch, read only right after a successful cycle, so it is not part of the snapshot.
    vector_t m_last_non_integral;
    pid_performance_tracker<Scalar, NY, Policies...> m_perf;
    bool m_first_step;
    bool m_integral_frozen;
    bool m_saturated;
    pid_health m_health;

    pid(unchecked_t, const config_type& cfg, const detail::pid_gains<Scalar, NY>& gains)
        : m_cfg{cfg}
        , m_gains{gains}
        , m_integral{vector_t::Zero()}
        , m_prev_error{vector_t::Zero()}
        , m_prev_derivative{vector_t::Zero()}
        , m_prev_meas{vector_t::Zero()}
        , m_prev_sp{vector_t::Zero()}
        , m_prev_output{vector_t::Zero()}
        , m_accumulated_output{vector_t::Zero()}
        , m_prev_ff{vector_t::Zero()}
        , m_filtered_sp{vector_t::Zero()}
        , m_filtered_meas{vector_t::Zero()}
        , m_prev_deriv_filtered{vector_t::Zero()}
        , m_constraint_direction{direction_t::Zero()}
        , m_last_non_integral{vector_t::Zero()}
        , m_first_step{true}
        , m_integral_frozen{false}
        , m_saturated{false}
        , m_health{pid_health::ok}
    {
    }

    cycle_state capture_cycle_state() const
    {
        return cycle_state{.integral = m_integral,
                .prev_error = m_prev_error,
                .prev_derivative = m_prev_derivative,
                .prev_meas = m_prev_meas,
                .prev_sp = m_prev_sp,
                .prev_output = m_prev_output,
                .accumulated_output = m_accumulated_output,
                .prev_ff = m_prev_ff,
                .filtered_sp = m_filtered_sp,
                .filtered_meas = m_filtered_meas,
                .prev_deriv_filtered = m_prev_deriv_filtered,
                .constraint_direction = m_constraint_direction,
                .performance = m_perf,
                .first_step = m_first_step,
                .saturated = m_saturated};
    }

    void restore_cycle_state(const cycle_state& previous)
    {
        m_integral = previous.integral;
        m_prev_error = previous.prev_error;
        m_prev_derivative = previous.prev_derivative;
        m_prev_meas = previous.prev_meas;
        m_prev_sp = previous.prev_sp;
        m_prev_output = previous.prev_output;
        m_accumulated_output = previous.accumulated_output;
        m_prev_ff = previous.prev_ff;
        m_filtered_sp = previous.filtered_sp;
        m_filtered_meas = previous.filtered_meas;
        m_prev_deriv_filtered = previous.prev_deriv_filtered;
        m_constraint_direction = previous.constraint_direction;
        m_perf = previous.performance;
        m_first_step = previous.first_step;
        m_saturated = previous.saturated;
    }

    bool track(const vector_t& target)
    {
        if constexpr(detail::contains_v<velocity_form, Policies...>)
            m_accumulated_output = target;
        else
        {
            auto next_integral = (target - m_last_non_integral).eval();
            if(!next_integral.allFinite())
                return false;
            m_integral = std::move(next_integral);
            // The actuator received the tracked command, not the automatic one.
            m_prev_output = target;
            m_constraint_direction = direction_t::Zero();
        }
        return true;
    }

    expected<void, pid_step_error> check_step(const vector_t& sp, const vector_t& meas, Scalar dt) const
    {
        if(!carried_state_finite())
            return unexpected(pid_step_error::non_finite_state);
        // dt <= 0 alone would admit a NaN step, since every comparison against a NaN is false.
        if(!std::isfinite(dt) || dt <= Scalar{0})
            return unexpected(pid_step_error::invalid_timestep);
        if(!sp.allFinite())
            return unexpected(pid_step_error::non_finite_setpoint);
        if(!meas.allFinite())
            return unexpected(pid_step_error::non_finite_measurement);
        return {};
    }

    // Each group is gated on the policy that makes the member live, so a value seeded into a member the
    // composed controller never reads cannot cause a rejection the output would not have suffered from.
    bool carried_state_finite() const
    {
        if(!m_prev_error.allFinite() || !m_prev_meas.allFinite() || !m_prev_sp.allFinite())
            return false;
        if constexpr(detail::contains_v<velocity_form, Policies...>)
        {
            if(!m_accumulated_output.allFinite() || !m_prev_derivative.allFinite())
                return false;
        }
        else if(!m_integral.allFinite())
            return false;
        return policy_state_finite();
    }

    bool policy_state_finite() const
    {
        return (!detail::contains_v<rate_limit, Policies...> || m_prev_output.allFinite()) && (!detail::has_policy_v<feed_forward, Policies...> || m_prev_ff.allFinite()) &&
                (!detail::contains_v<setpoint_filter, Policies...> || m_filtered_sp.allFinite()) && (!detail::contains_v<pv_filter, Policies...> || m_filtered_meas.allFinite()) &&
                (!detail::contains_v<deriv_filter, Policies...> || m_prev_deriv_filtered.allFinite());
    }

    pid_step_error latch_health(pid_step_error fault)
    {
        if(fault == pid_step_error::non_finite_state)
            m_health = pid_health::non_finite_carried_state;
        return fault;
    }

    vector_t apply_setpoint_filter(const vector_t& sp, Scalar dt)
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

    vector_t apply_pv_filter(const vector_t& meas, Scalar dt)
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

    vector_t compute_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt)
    {
        if constexpr(detail::contains_v<velocity_form, Policies...>)
            return compute_velocity_form(e, sp, filtered_sp, filtered_meas, dt);
        else
            return compute_position_form(e, sp, filtered_sp, filtered_meas, dt);
    }

    // The incremental form as the first difference of the positional terms, so step 1 has no
    // derivative increment and each cycle's own dt divides its own difference.
    vector_t compute_velocity_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt)
    {
        auto dp = (compute_proportional_term(filtered_sp, filtered_meas) - compute_proportional_term(m_prev_sp, m_prev_meas)).eval();
        auto d = compute_derivative_term(filtered_sp, filtered_meas, dt);
        auto delta_u = (dp + integral_increment(e, dt) + (d - m_prev_derivative)).eval();
        m_prev_derivative = d;
        auto emitted = accumulate_output(clip_increment(apply_feed_forward_velocity(delta_u, sp, dt), dt).value);
        update_state(e, filtered_meas, filtered_sp, emitted);
        return emitted;
    }

    // Clamping the accumulated output rather than the increment keeps it in range while still letting it
    // rise and fall against an asymmetric limit. The clamp returns its input when inside, so != is exact.
    vector_t accumulate_output(const vector_t& delta_u)
    {
        auto target = (m_accumulated_output + delta_u).eval();
        auto clamped = target.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        m_saturated = (clamped.array() != target.array()).any();
        auto emitted = (clamped - m_accumulated_output).eval();
        m_accumulated_output = clamped;
        return emitted;
    }

    // The velocity form emits increments, so it injects the change in the feed-forward level, not the level:
    // a constant feed-forward then contributes nothing to the increment and the actuator does not drift.
    vector_t apply_feed_forward_velocity(vector_t delta_u, const vector_t& sp, Scalar dt)
    {
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                auto ff = m_cfg.template policy<ff_policy_t>().feed_forward_cb(sp, dt).eval();
                auto delta_ff = (ff - m_prev_ff).eval();
                m_prev_ff = ff;
                return (delta_u + delta_ff).eval();
            }
        }
        return delta_u;
    }

    // @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006, Ch. 3-4
    vector_t compute_position_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt)
    {
        auto p = compute_proportional_term(filtered_sp, filtered_meas);
        auto integral = compute_integral_term(e, dt);
        auto d = compute_derivative_term(filtered_sp, filtered_meas, dt);
        auto u_unconstrained = compute_raw_output(p, integral, d, sp, dt);
        auto limited = apply_rate_limit(u_unconstrained, dt);
        auto u_sat = limited.value.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        auto saturated = (limited.value.array() < m_cfg.output_min.array() || limited.value.array() > m_cfg.output_max.array()).eval();
        m_saturated = saturated.any();
        // Anti-windup feeds back against the command taken before the rate limiter, so a rate-limited ramp,
        // where the applied output already equals the limited command, engages it as well as saturation does.
        apply_anti_windup(u_sat, u_unconstrained, (limited.clipped || saturated).eval(), dt);
        update_state(e, filtered_meas, filtered_sp, u_sat);
        return u_sat;
    }

    // @cite astrom2006 -- Ch. 3.5 (setpoint weighting b parameter)
    vector_t compute_proportional_term(const vector_t& filtered_sp, const vector_t& filtered_meas)
    {
        vector_t ep;
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            ep[i] = m_cfg.b[i] * filtered_sp[i] - filtered_meas[i];
        return m_gains.kp.cwiseProduct(ep).eval();
    }

    // @cite astrom2006 -- Ch. 3.3 (integral action, discretization methods)
    vector_t integral_increment(const vector_t& e, Scalar dt) const
    {
        if(m_integral_frozen)
            return vector_t::Zero();
        if constexpr(detail::contains_v<forward_euler, Policies...>)
            return (m_gains.ki.cwiseProduct(m_prev_error) * dt).eval();
        else if constexpr(detail::contains_v<tustin, Policies...>)
            return (m_gains.ki.cwiseProduct(((e + m_prev_error) * Scalar{0.5}).eval()) * dt).eval();
        else
            return (m_gains.ki.cwiseProduct(e) * dt).eval();
    }

    vector_t compute_integral_term(const vector_t& e, Scalar dt)
    {
        if(!m_integral_frozen)
            m_integral = (m_integral + gate_increment(integral_increment(e, dt), e)).eval();
        return m_integral;
    }

    // Clamping and conditional integration decide before the output is formed, so no increment that formed a
    // command is removed afterward; clamping therefore reads the previous cycle's constraint.
    vector_t gate_increment(vector_t increment, const vector_t& e) const
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if constexpr(detail::contains_v<anti_windup<clamping>, Policies...>)
            {
                const std::int8_t direction = m_constraint_direction[i];
                if((direction > 0 && increment[i] > Scalar{0}) || (direction < 0 && increment[i] < Scalar{0}))
                    increment[i] = Scalar{0};
            }
            else if constexpr(detail::contains_v<anti_windup<conditional_integration>, Policies...>)
            {
                const auto& threshold = m_cfg.template policy<anti_windup<conditional_integration>>().error_threshold;
                if(std::abs(e[i]) > threshold[static_cast<std::size_t>(i)])
                    increment[i] = Scalar{0};
            }
        }
        return increment;
    }

    // @cite astrom2006 -- Ch. 3.4 (derivative action, setpoint weighting c parameter)
    vector_t compute_derivative_term(const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt)
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

    vector_t compute_raw_derivative(const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt)
    {
        if(m_cfg.derivative_on_error)
        {
            vector_t ed_curr, ed_prev;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                ed_curr[i] = m_cfg.c[i] * filtered_sp[i] - filtered_meas[i];
                ed_prev[i] = m_cfg.c[i] * m_prev_sp[i] - m_prev_meas[i];
            }
            return m_gains.kd.cwiseProduct((ed_curr - ed_prev) / dt).eval();
        }
        else
        {
            auto dm_dt = ((filtered_meas - m_prev_meas) / dt).eval();
            return (-m_gains.kd.cwiseProduct(dm_dt)).eval();
        }
    }

    vector_t apply_derivative_filter(vector_t d, Scalar dt)
    {
        if constexpr(detail::contains_v<deriv_filter, Policies...>)
        {
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                if(m_gains.tf[i] > Scalar{0})
                {
                    Scalar alpha = m_gains.tf[i] / (m_gains.tf[i] + dt);
                    d[i] = alpha * m_prev_deriv_filtered[i] + (Scalar{1} - alpha) * d[i];
                }
            m_prev_deriv_filtered = d;
        }
        return d;
    }

    vector_t compute_raw_output(const vector_t& p, const vector_t& integral, const vector_t& d, const vector_t& sp, Scalar dt)
    {
        auto u_raw = (p + integral + d).eval();
        m_last_non_integral = (p + d).eval();
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                // The feed-forward callable may be stateful, so it is evaluated once per cycle.
                const vector_t ff = m_cfg.template policy<ff_policy_t>().feed_forward_cb(sp, dt);
                m_last_non_integral += ff;
                return (u_raw + ff).eval();
            }
        }
        return u_raw;
    }

    limited_command apply_rate_limit(const vector_t& u_raw, Scalar dt) const
    {
        if constexpr(!detail::contains_v<rate_limit, Policies...>)
            return {u_raw, mask_t::Constant(false)};
        auto limited = clip_increment((u_raw - m_prev_output).eval(), dt);
        // Re-forming p + (u - p) rounds, so a channel within its bound keeps its input bit for bit.
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            limited.value[i] = limited.clipped[i] ? m_prev_output[i] + limited.value[i] : u_raw[i];
        return limited;
    }

    // A zero rate leaves its channel unlimited, and an infinite one gives an infinite bound that never binds.
    limited_command clip_increment(const vector_t& delta, Scalar dt) const
    {
        limited_command limited{delta, mask_t::Constant(false)};
        if constexpr(detail::contains_v<rate_limit, Policies...>)
        {
            const auto& rate_max = m_cfg.template policy<rate_limit>().rate_max;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                const Scalar rate = rate_max[static_cast<std::size_t>(i)];
                const Scalar step = rate * dt;
                if(rate > Scalar{0} && (delta[i] > step || delta[i] < -step))
                {
                    limited.value[i] = delta[i] > step ? step : -step;
                    limited.clipped[i] = true;
                }
            }
        }
        return limited;
    }

    // A frozen integrator is held against back-calculation too, while the clamping verdict is recorded regardless,
    // since it describes the emitted command rather than the integrator.
    void apply_anti_windup(const vector_t& u_sat, const vector_t& u_unconstrained, const mask_t& constrained, Scalar dt)
    {
        if constexpr(detail::contains_v<anti_windup<clamping>, Policies...>)
            m_constraint_direction = constraint_direction(u_sat, u_unconstrained, constrained);
        else if constexpr(detail::contains_v<anti_windup<back_calc>, Policies...>)
        {
            if(!m_integral_frozen)
                apply_back_calculation(u_sat, u_unconstrained, dt);
        }
    }

    void apply_back_calculation(const vector_t& u_sat, const vector_t& u_unconstrained, Scalar dt)
    {
        auto feedback = m_gains.kb.cwiseProduct((u_sat - u_unconstrained).eval()).eval();
        m_integral = (m_integral + feedback * dt).eval();
    }

    // @cite caparroz2026 -- clamping decided from the previous sample's saturation error, Algorithm 2. The
    // verdict is taken from the masks; the unconstrained command only names the side of the violated limit.
    direction_t constraint_direction(const vector_t& u_sat, const vector_t& u_unconstrained, const mask_t& constrained) const
    {
        direction_t direction = direction_t::Zero();
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            if(constrained[i] && u_unconstrained[i] != u_sat[i])
                direction[i] = u_unconstrained[i] > u_sat[i] ? std::int8_t{1} : std::int8_t{-1};
        return direction;
    }

    void update_state(const vector_t& e, const vector_t& filtered_meas, const vector_t& filtered_sp, const vector_t& output)
    {
        m_prev_error = e;
        m_prev_meas = filtered_meas;
        m_prev_sp = filtered_sp;
        m_prev_output = output;
        m_first_step = false;
    }
};

}

#endif
