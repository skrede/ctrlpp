#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

/// @brief Policy-based PID controller with compile-time feature composition.
///
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/pid_error.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"
#include "ctrlpp/control/pid_performance.h"

#include "ctrlpp/control/detail/pid_gains.h"

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

public:
    using config_type = pid_config<Scalar, NY, Policies...>;
    using vector_t = Vector<Scalar, NY>;

    /// @brief Construct a controller, or report why the configuration cannot work.
    /// This is the only construction path.
    ///
    /// Rejections, checked in order:
    ///  * a NaN or infinite kp, ki or kd component, or in ISA form a NaN or
    ///    negative-infinite integral time -> pid_config_error::non_finite_gain
    ///  * in ISA form a zero integral time -> pid_config_error::zero_integral_time
    ///  * a NaN output limit, output_min above output_max, output_min at
    ///    +infinity or output_max at -infinity -> pid_config_error::invalid_output_limits
    ///  * a NaN or infinite b or c component -> pid_config_error::non_finite_setpoint_weight
    ///  * a setpoint_filter or pv_filter time constant that is NaN, negative or
    ///    +infinity -> pid_config_error::invalid_filter_time_constant
    ///  * a NaN or negative rate_limit rate_max -> pid_config_error::invalid_rate_limit
    ///  * a NaN or infinite explicit back-calculation kb ->
    ///    pid_config_error::invalid_back_calculation_gain
    ///  * a NaN conditional-integration error_threshold -> pid_config_error::invalid_error_threshold
    ///  * a NaN or infinite deriv_filter n, or an explicit deriv_filter tf that is
    ///    not finite and positive -> pid_config_error::invalid_derivative_filter
    ///  * with oscillation_detect, an amplitude or ultimate_frequency component that
    ///    is not positive and finite -> pid_config_error::missing_oscillation_parameters
    ///  * with oscillation_detect, a load_count_limit or supervision_periods that is
    ///    not positive and finite -> pid_config_error::invalid_oscillation_supervision
    ///  * a default back-calculation kb on a channel with integral action whose
    ///    nonzero gains differ in sign -> pid_config_error::undefined_tracking_time
    ///  * a derivative channel without an explicit tf whose kp is zero, whose n is
    ///    not positive, or whose kd and kp differ in sign
    ///    -> pid_config_error::undefined_derivative_filter
    ///  * an ISA-form internal gain, kp over the integral time or kp times the
    ///    derivative time, that overflows, a default filter time constant
    ///    kd/(kp*n) that is not positive and finite, or a default back-calculation
    ///    kb that is not finite -> pid_config_error::unrepresentable_gain
    static auto create(const config_type& cfg) -> expected<pid, pid_config_error>
    {
        auto derived = detail::derive_pid_gains(cfg);
        if(!derived)
            return unexpected(derived.error());
        return pid{unchecked_t{}, cfg, *derived};
    }

    /// @brief Produce the control command for one cycle.
    ///
    /// The cycle is classified before any member is written, so a rejected
    /// cycle leaves every piece of carried state -- the integrator, the error
    /// and derivative histories, the input filter states, the accumulated output -- bitwise
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
    /// Configuration is checked where it is given, by `create`, `set_params`
    /// and `set_integral`, and a rejected cycle mutates nothing.
    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt) -> expected<vector_t, pid_step_error>
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

    /// @brief Produce the control command for one cycle and re-seed the
    /// controller so the output tracks an external signal (bumpless transfer).
    ///
    /// The returned value is the cycle's own computed command, or increment in
    /// the velocity form. In the position form the integrator becomes the
    /// tracking signal clamped to the output limits, minus the cycle's
    /// proportional, derivative and feed-forward terms taken before rate
    /// limiting and saturation, so it never holds an offset that drives the
    /// unconstrained command past a limit. In the velocity form the accumulated
    /// output becomes the clamped tracking signal.
    ///
    /// In the position form the rate limiter's anchor also becomes the clamped
    /// tracking signal and the carried clamping verdict is cleared, so the next
    /// automatic cycle neither rate-limits from the discarded automatic command
    /// nor withholds its increment on that command's constraint.
    ///
    /// Tracking assigns the integrator even while `freeze_integral` holds it:
    /// like `set_integral` it is an explicit assignment, and the freeze holds the
    /// integrator only against the integral increment and anti-windup.
    ///
    /// The tracking signal is checked BEFORE the cycle is delegated, because it
    /// is written straight into carried state once the cycle succeeds:
    /// admitting a non-finite one would destroy that state after the cycle
    /// had already committed its other state, which is exactly the
    /// partially-applied step the reject-before-mutate rule exists to forbid.
    ///
    /// A rejection from the delegated cycle is forwarded unchanged, carrying
    /// its own specific cause, and the tracking assignment is not performed --
    /// so a rejected cycle leaves the carried state bitwise unchanged here too.
    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt, const vector_t& tracking_signal) -> expected<vector_t, pid_step_error>
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

    /// @brief Replace the configuration, or report why it cannot work and keep
    /// the current one bitwise. Refusals are those of `create`, in its order.
    auto set_params(const config_type& new_cfg) -> expected<void, pid_config_error>
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
        m_prev_derivative = vector_t::Zero();
        m_filtered_sp = vector_t::Zero();
        m_filtered_meas = vector_t::Zero();
        m_prev_deriv_filtered = vector_t::Zero();
        m_first_step = true;
        m_integral_frozen = false;
        m_saturated = false;
        m_constraint_direction = direction_t::Zero();
        // Every member the latched status describes has just been replaced, so
        // the status no longer describes anything and clearing it is honest.
        m_health = pid_health::ok;
        if constexpr(detail::has_policy_v<perf_assessment, Policies...>)
            m_perf.reset();
    }

    void freeze_integral(bool freeze = true) { m_integral_frozen = freeze; }

    auto set_integral(const vector_t& val) -> expected<void, pid_config_error>
    {
        if(!val.allFinite())
            return unexpected(pid_config_error::non_finite_integral);
        m_integral = val;
        return {};
    }

    template <typename Metric>
    auto metric() const -> const vector_t&
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_perf.template metric<Metric>();
    }

    auto oscillating() const -> bool
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
    pid(unchecked_t, const config_type& cfg, const detail::pid_gains<Scalar, NY>& gains) : m_cfg{cfg}, m_gains{gains}, m_constraint_direction{direction_t::Zero()},
        m_last_non_integral{vector_t::Zero()}
    {
    }

    struct cycle_state
    {
        vector_t integral;
        vector_t prev_error;
        vector_t prev_derivative;
        vector_t prev_meas;
        vector_t prev_sp;
        vector_t prev_output;
        vector_t accumulated_output;
        vector_t prev_ff;
        vector_t filtered_sp;
        vector_t filtered_meas;
        vector_t prev_deriv_filtered;
        direction_t constraint_direction;
        pid_performance_tracker<Scalar, NY, Policies...> performance;
        bool first_step;
        bool saturated;
    };

    auto capture_cycle_state() const -> cycle_state
    {
        return cycle_state{
            .integral = m_integral,
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

    auto track(const vector_t& target) -> bool
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
            if(!m_accumulated_output.allFinite() || !m_prev_derivative.allFinite())
                return false;
        }
        else if(!m_integral.allFinite())
            return false;
        return policy_state_finite();
    }

    auto policy_state_finite() const -> bool
    {
        return (!detail::contains_v<rate_limit, Policies...> || m_prev_output.allFinite())
            && (!detail::has_policy_v<feed_forward, Policies...> || m_prev_ff.allFinite())
            && (!detail::contains_v<setpoint_filter, Policies...> || m_filtered_sp.allFinite())
            && (!detail::contains_v<pv_filter, Policies...> || m_filtered_meas.allFinite())
            && (!detail::contains_v<deriv_filter, Policies...> || m_prev_deriv_filtered.allFinite());
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

    auto compute_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
    {
        if constexpr(detail::contains_v<velocity_form, Policies...>)
            return compute_velocity_form(e, sp, filtered_sp, filtered_meas, dt);
        else
            return compute_position_form(e, sp, filtered_sp, filtered_meas, dt);
    }

    // The incremental form as the first difference of the positional terms, so step 1 has no
    // derivative increment and each cycle's own dt divides its own difference.
    auto compute_velocity_form(const vector_t& e, const vector_t& sp, const vector_t& filtered_sp, const vector_t& filtered_meas, Scalar dt) -> vector_t
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
    auto accumulate_output(const vector_t& delta_u) -> vector_t
    {
        auto target = (m_accumulated_output + delta_u).eval();
        auto clamped = target.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        m_saturated = (clamped.array() != target.array()).any();
        auto emitted = (clamped - m_accumulated_output).eval();
        m_accumulated_output = clamped;
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
        auto integral = compute_integral_term(e, dt);
        auto d = compute_derivative_term(filtered_sp, filtered_meas, dt);
        // Unconstrained control command, before rate limiting and output saturation.
        auto u_unconstrained = compute_raw_output(p, integral, d, sp, dt);
        auto limited = apply_rate_limit(u_unconstrained, dt);
        auto u_sat = limited.value.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval();
        auto saturated = (limited.value.array() < m_cfg.output_min.array() || limited.value.array() > m_cfg.output_max.array()).eval();
        m_saturated = saturated.any();
        // Anti-windup feeds back against the command taken before the rate limiter, so a
        // rate-limited ramp, where the applied output already equals the limited command,
        // engages it as well as output saturation does.
        apply_anti_windup(u_sat, u_unconstrained, (limited.clipped || saturated).eval(), dt);
        update_state(e, filtered_meas, filtered_sp, u_sat);
        return u_sat;
    }

    /// @cite astrom2006 -- Ch. 3.5 (setpoint weighting b parameter)
    auto compute_proportional_term(const vector_t& filtered_sp, const vector_t& filtered_meas) -> vector_t
    {
        vector_t ep;
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            ep[i] = m_cfg.b[i] * filtered_sp[i] - filtered_meas[i];
        return m_gains.kp.cwiseProduct(ep).eval();
    }

    /// @cite astrom2006 -- Ch. 3.3 (integral action, discretization methods)
    auto integral_increment(const vector_t& e, Scalar dt) const -> vector_t
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

    auto compute_integral_term(const vector_t& e, Scalar dt) -> vector_t
    {
        if(!m_integral_frozen)
            m_integral = (m_integral + gate_increment(integral_increment(e, dt), e)).eval();
        return m_integral;
    }

    // Clamping and conditional integration decide before the output is formed, so no increment that formed a
    // command is removed afterward; clamping therefore reads the previous cycle's constraint.
    auto gate_increment(vector_t increment, const vector_t& e) const -> vector_t
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
            return m_gains.kd.cwiseProduct((ed_curr - ed_prev) / dt).eval();
        }
        else
        {
            auto dm_dt = ((filtered_meas - m_prev_meas) / dt).eval();
            return (-m_gains.kd.cwiseProduct(dm_dt)).eval();
        }
    }

    auto apply_derivative_filter(vector_t d, Scalar dt) -> vector_t
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

    auto compute_raw_output(const vector_t& p, const vector_t& integral, const vector_t& d, const vector_t& sp, Scalar dt) -> vector_t
    {
        auto u_raw = (p + integral + d).eval();
        m_last_non_integral = (p + d).eval();
        if constexpr(detail::has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                // The feed-forward callable may be stateful, so it is evaluated once per cycle.
                const vector_t ff = m_cfg.template policy<ff_policy_t>().ff_func(sp, dt);
                m_last_non_integral += ff;
                return (u_raw + ff).eval();
            }
        }
        return u_raw;
    }

    auto apply_rate_limit(const vector_t& u_raw, Scalar dt) const -> limited_command
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
    auto clip_increment(const vector_t& delta, Scalar dt) const -> limited_command
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

    /// @cite caparroz2026 -- clamping decided from the previous sample's saturation error, Algorithm 2. The
    /// verdict is taken from the masks; the unconstrained command only names the side of the violated limit.
    auto constraint_direction(const vector_t& u_sat, const vector_t& u_unconstrained, const mask_t& constrained) const -> direction_t
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

    config_type m_cfg;
    detail::pid_gains<Scalar, NY> m_gains;
    vector_t m_integral = vector_t::Zero();
    vector_t m_prev_error = vector_t::Zero();
    vector_t m_prev_derivative = vector_t::Zero();
    vector_t m_prev_meas = vector_t::Zero();
    vector_t m_prev_sp = vector_t::Zero();
    vector_t m_prev_output = vector_t::Zero();
    vector_t m_accumulated_output = vector_t::Zero();
    vector_t m_prev_ff = vector_t::Zero();
    vector_t m_filtered_sp = vector_t::Zero();
    vector_t m_filtered_meas = vector_t::Zero();
    vector_t m_prev_deriv_filtered = vector_t::Zero();
    direction_t m_constraint_direction;
    // Per-cycle scratch, read only right after a successful cycle, so it is not part of the snapshot.
    vector_t m_last_non_integral;
    pid_performance_tracker<Scalar, NY, Policies...> m_perf;
    bool m_first_step{true};
    bool m_integral_frozen{false};
    bool m_saturated{false};
    pid_health m_health{pid_health::ok};
};

}

#endif
