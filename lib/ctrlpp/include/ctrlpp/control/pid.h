#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

/// @brief Policy-based PID controller with compile-time feature composition.
///
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"
#include "ctrlpp/control/pid_performance.h"

#include <cmath>
#include <limits>

namespace ctrlpp
{

/// @brief Structured failure modes of a single `pid::compute` cycle.
///
/// Every enumerator is an exact domain condition, not a tuning preference. The
/// order below is the order the guard tests them, and it is a severity order:
/// the controller's own carried state first, then the caller's clock, then the
/// two signals. The rule is that the most upstream cause is named, because a
/// caller told "your measurement is bad" would replace a working sensor while
/// the real fault sits in the integrator that a previous cycle poisoned.
///
///  * non_finite_state       : some piece of carried state -- the integrator,
///                             the error history, an input filter state, the
///                             accumulated output -- is already non-finite when
///                             the cycle begins. Nothing this cycle can produce
///                             is meaningful, whatever the arguments are. Only
///                             `reset` (or `set_integral` for the integrator
///                             alone) recovers.
///  * invalid_timestep       : the step is not a positive finite duration. A
///                             non-positive step is a clock that did not advance
///                             or ran backwards; a non-finite one is a clock
///                             that produced garbage. Both are the caller's
///                             timing, not the caller's signals, and both are
///                             fatal to the cycle: the step divides the
///                             derivative and multiplies the integral, so it
///                             poisons the command even from a perfect setpoint
///                             and measurement.
///  * non_finite_setpoint    : the commanded setpoint has a non-finite
///                             component. The repair is in whatever generates
///                             the reference -- a trajectory, an outer loop, an
///                             operator input.
///  * non_finite_measurement : the process variable has a non-finite component.
///                             The repair is in the sensor or the estimator
///                             feeding it. This is a genuinely different fault
///                             from a bad setpoint even though the error is
///                             their difference: the two arrive from different
///                             subsystems and are fixed in different places.
///  * non_finite_tracking_signal : the external tracking signal handed to the
///                             four-argument overload has a non-finite
///                             component. It is checked before the cycle runs
///                             because it is back-assigned into the integrator,
///                             so admitting it would destroy the integrator
///                             after an otherwise valid cycle had already
///                             committed.
enum class pid_step_error
{
    non_finite_state,
    invalid_timestep,
    non_finite_setpoint,
    non_finite_measurement,
    non_finite_tracking_signal,
    non_finite_result,
};

/// @brief Persistent state-health status of a `pid`.
///
/// A per-cycle result cannot answer whether the controller is still carrying
/// damage from a cycle several samples ago, because that question outlives the
/// call. The status latches until `reset` replaces the state it describes.
///
///  * ok                      : every cycle so far began from finite carried
///                              state.
///  * non_finite_carried_state : a cycle found the carried state already
///                              non-finite. The usual route in is a finite but
///                              extreme configuration -- an infinite gain
///                              multiplies a finite error into an infinite
///                              command, which `update_state` then stores --
///                              or a non-finite value seeded through
///                              `set_integral` or `set_params`, neither of
///                              which is fallible. A rejected cycle does NOT
///                              set this: a rejection mutates nothing, so it
///                              leaves the controller exactly as healthy as it
///                              was.
enum class pid_health
{
    ok,
    non_finite_carried_state,
};

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

    /// @brief Produce the control command for one cycle.
    ///
    /// The cycle is classified before any member is written, so a rejected
    /// cycle leaves every piece of carried state -- the integrator, both error
    /// histories, the input filter states, the accumulated output -- bitwise
    /// unchanged, and the caller may retry on the next sample.
    ///
    /// **What a refusal means for the caller.** A rejected cycle produced NO
    /// command. That is a materially different situation from an estimator
    /// declining to fold in a measurement, where the estimate simply stands: an
    /// actuator is going to be driven by something regardless of what this
    /// function returns. The caller must choose that something, and this
    /// controller deliberately does not choose for it, because the right choice
    /// is a property of the plant and not of the controller: holding the last
    /// command is correct for a slow thermal loop and dangerous for an unstable
    /// attitude loop, where zero or a fail-over path is correct instead. The
    /// three defensible responses are to hold the command the last successful
    /// cycle returned, to command a configured safe value, or to fail over to a
    /// redundant channel.
    ///
    /// What the controller will no longer do is hand back the previous command
    /// dressed as a fresh one. On a non-positive step it used to return the
    /// stored output on the success path, which a caller could not distinguish
    /// from a command the controller had actually computed -- so a stopped
    /// clock read as a steady loop.
    ///
    /// `set_params`, `set_integral`, `freeze_integral` and `reset` are NOT
    /// fallible and are not converted here. A non-finite gain or a non-finite
    /// seeded integrator therefore still enters through them; the first cycle
    /// that follows rejects with `non_finite_state` and latches `health()`,
    /// which is what that query exists for.
    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt) -> expected<vector_t, pid_step_error>
    {
        if(const auto step = check_step(sp, meas, dt); !step)
            return unexpected(latch_health(step.error()));

        auto const previous = capture_cycle_state();
        auto filtered_sp = apply_setpoint_filter(sp, dt);
        auto filtered_meas = apply_pv_filter(meas, dt);
        auto e = (filtered_sp - filtered_meas).eval();
        m_perf.accumulate(e, dt);

        auto result = [&] {
            if constexpr(detail::contains_v<velocity_form, Policies...>)
                return compute_velocity_form(
                    e, sp, filtered_sp, filtered_meas, dt);
            else
                return compute_position_form(
                    e, sp, filtered_sp, filtered_meas, dt);
        }();

        if(!result.allFinite() || !carried_state_finite()
            || !m_perf.all_finite())
        {
            restore_cycle_state(previous);
            return unexpected(pid_step_error::non_finite_result);
        }
        return result;
    }

    /// @brief Produce the control command for one cycle and back-assign the
    /// integrator so the output tracks an external signal (bumpless transfer).
    ///
    /// The tracking signal is checked BEFORE the cycle is delegated, because it
    /// is written straight into the integrator once the cycle succeeds:
    /// admitting a non-finite one would destroy the integrator after the cycle
    /// had already committed its other state, which is exactly the
    /// partially-applied step the reject-before-mutate rule exists to forbid.
    ///
    /// A rejection from the delegated cycle is forwarded unchanged, carrying
    /// its own specific cause, and the tracking assignment is not performed --
    /// so a rejected cycle leaves the integrator bitwise unchanged here too.
    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt, const vector_t& tracking_signal) -> expected<vector_t, pid_step_error>
    {
        if(!tracking_signal.allFinite())
            return unexpected(pid_step_error::non_finite_tracking_signal);

        auto const previous = capture_cycle_state();
        auto u = compute(sp, meas, dt);
        if(!u)
            return u;
        if constexpr(!detail::contains_v<velocity_form, Policies...>)
        {
            auto non_integral = (*u - m_integral).eval();
            auto next_integral = (tracking_signal - non_integral).eval();
            if(!next_integral.allFinite())
            {
                restore_cycle_state(previous);
                return unexpected(pid_step_error::non_finite_result);
            }
            m_integral = std::move(next_integral);
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

    /// @brief Report whether the controller is still carrying non-finite state
    /// from an earlier cycle. Latches until `reset`; a rejected cycle does not
    /// set it.
    auto health() const -> pid_health { return m_health; }

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
        // Every member the latched status describes has just been replaced, so
        // the status no longer describes anything and clearing it is honest.
        m_health = pid_health::ok;
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
    struct cycle_state
    {
        vector_t integral;
        vector_t prev_error;
        vector_t prev_prev_error;
        vector_t prev_meas;
        vector_t prev_sp;
        vector_t prev_output;
        vector_t accumulated_output;
        vector_t prev_ff;
        vector_t filtered_sp;
        vector_t filtered_meas;
        vector_t prev_deriv_filtered;
        pid_performance_tracker<Scalar, NY, Policies...> performance;
        bool first_step;
        bool saturated;
    };

    auto capture_cycle_state() const -> cycle_state
    {
        return cycle_state{
            .integral = m_integral,
            .prev_error = m_prev_error,
            .prev_prev_error = m_prev_prev_error,
            .prev_meas = m_prev_meas,
            .prev_sp = m_prev_sp,
            .prev_output = m_prev_output,
            .accumulated_output = m_accumulated_output,
            .prev_ff = m_prev_ff,
            .filtered_sp = m_filtered_sp,
            .filtered_meas = m_filtered_meas,
            .prev_deriv_filtered = m_prev_deriv_filtered,
            .performance = m_perf,
            .first_step = m_first_step,
            .saturated = m_saturated};
    }

    void restore_cycle_state(const cycle_state& previous)
    {
        m_integral = previous.integral;
        m_prev_error = previous.prev_error;
        m_prev_prev_error = previous.prev_prev_error;
        m_prev_meas = previous.prev_meas;
        m_prev_sp = previous.prev_sp;
        m_prev_output = previous.prev_output;
        m_accumulated_output = previous.accumulated_output;
        m_prev_ff = previous.prev_ff;
        m_filtered_sp = previous.filtered_sp;
        m_filtered_meas = previous.filtered_meas;
        m_prev_deriv_filtered = previous.prev_deriv_filtered;
        m_perf = previous.performance;
        m_first_step = previous.first_step;
        m_saturated = previous.saturated;
    }

    /// @brief Classify a cycle's operands without touching a single member.
    ///
    /// The order is the severity order documented on `pid_step_error`: the
    /// carried state first, then the step, then the two signals. The cost is
    /// one finiteness scan of each live operand -- a fixed number of NY-element
    /// reads plus one scalar test, no data-dependent branching and no
    /// allocation, NY being a compile-time template parameter.
    auto check_step(const vector_t& sp, const vector_t& meas, Scalar dt) const -> expected<void, pid_step_error>
    {
        if(!carried_state_finite())
            return unexpected(pid_step_error::non_finite_state);
        // isfinite rejects both a NaN step and an infinite one; the second test
        // then rejects zero and negative. A NaN cannot be caught by dt <= 0
        // alone, because every comparison against a NaN is false -- which is
        // precisely how a NaN step used to walk past the old guard.
        if(!std::isfinite(dt) || dt <= Scalar{0})
            return unexpected(pid_step_error::invalid_timestep);
        if(!sp.allFinite())
            return unexpected(pid_step_error::non_finite_setpoint);
        if(!meas.allFinite())
            return unexpected(pid_step_error::non_finite_measurement);
        return {};
    }

    /// @brief Test every member that can reach this cycle's command.
    ///
    /// Each group is gated on the policy that makes the member live, so a value
    /// seeded into a member the composed controller never reads cannot cause a
    /// rejection the output would not have suffered from.
    auto carried_state_finite() const -> bool
    {
        if(!m_prev_error.allFinite() || !m_prev_meas.allFinite() || !m_prev_sp.allFinite())
            return false;
        if constexpr(detail::contains_v<velocity_form, Policies...>)
        {
            if(!m_accumulated_output.allFinite() || !m_prev_prev_error.allFinite())
                return false;
        }
        else
        {
            if(!m_integral.allFinite())
                return false;
        }
        if constexpr(detail::contains_v<rate_limit, Policies...>)
        {
            if(!m_prev_output.allFinite())
                return false;
        }
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            if(!m_prev_ff.allFinite())
                return false;
        }
        if constexpr(detail::contains_v<setpoint_filter, Policies...>)
        {
            if(!m_filtered_sp.allFinite())
                return false;
        }
        if constexpr(detail::contains_v<pv_filter, Policies...>)
        {
            if(!m_filtered_meas.allFinite())
                return false;
        }
        if constexpr(detail::contains_v<deriv_filter, Policies...>)
        {
            if(!m_prev_deriv_filtered.allFinite())
                return false;
        }
        return true;
    }

    /// @brief Latch the persistent status for the one fault that describes the
    /// controller rather than the cycle's arguments, and pass the fault through
    /// so the caller still receives the specific cause on the failure channel.
    auto latch_health(pid_step_error fault) -> pid_step_error
    {
        if(fault == pid_step_error::non_finite_state)
            m_health = pid_health::non_finite_carried_state;
        return fault;
    }

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

    /// @cite astrom2006 -- Ch. 3.3 (integral action, discretization methods)
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
    pid_health m_health{pid_health::ok};
};

}

#endif
