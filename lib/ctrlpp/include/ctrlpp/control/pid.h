#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_H

// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/control/pid_error.h"
#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"
#include "ctrlpp/control/pid_performance.h"

#include "ctrlpp/control/detail/pid_gains.h"
#include "ctrlpp/control/detail/pid_cycle_state.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <cmath>
#include <cstddef>

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

        auto const previous = m_state;
        auto result = m_state.advance(m_cfg, m_gains, m_integral_frozen, sp, meas, dt, m_last_non_integral);
        if(!result.allFinite() || !m_state.finite() || !m_state.performance().all_finite())
        {
            m_state = previous;
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

        auto const previous = m_state;
        auto u = compute(sp, meas, dt);
        if(!u)
            return u;
        if(!m_state.track(tracking_signal.cwiseMax(m_cfg.output_min).cwiseMin(m_cfg.output_max).eval(), m_last_non_integral))
        {
            m_state = previous;
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

    const vector_t& error() const { return m_state.error(); }
    const vector_t& integral() const { return m_state.integral(); }
    const config_type& params() const { return m_cfg; }
    bool saturated() const { return m_state.saturated(); }

    pid_health health() const { return m_health; }

    void reset()
    {
        m_state = detail::pid_cycle_state<Scalar, NY, Policies...>{};
        m_integral_frozen = false;
        m_health = pid_health::ok;
    }

    void freeze_integral(bool freeze = true) { m_integral_frozen = freeze; }

    expected<void, pid_config_error> set_integral(const vector_t& val)
        requires(!detail::contains_v<velocity_form, Policies...>)
    {
        if(!val.allFinite())
            return unexpected(pid_config_error::non_finite_integral);
        m_state.assign_integral(val);
        return {};
    }

    template <typename Metric>
    const vector_t& metric() const
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_state.performance().template metric<Metric>();
    }

    bool oscillating() const
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        return m_state.performance().oscillating(m_cfg);
    }

    void reset_metrics()
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        m_state.performance().reset();
    }

private:
    config_type m_cfg;
    pid_health m_health;
    bool m_integral_frozen;
    // Per-cycle scratch, read only right after a successful cycle, so it is not part of the snapshot.
    vector_t m_last_non_integral;
    detail::pid_gains<Scalar, NY> m_gains;
    detail::pid_cycle_state<Scalar, NY, Policies...> m_state;

    pid(unchecked_t, const config_type& cfg, const detail::pid_gains<Scalar, NY>& gains)
        : m_cfg{cfg}
        , m_health{pid_health::ok}
        , m_integral_frozen{false}
        , m_last_non_integral{vector_t::Zero()}
        , m_gains{gains}
    {
    }

    expected<void, pid_step_error> check_step(const vector_t& sp, const vector_t& meas, Scalar dt) const
    {
        if(!m_state.finite())
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

    pid_step_error latch_health(pid_step_error fault)
    {
        if(fault == pid_step_error::non_finite_state)
            m_health = pid_health::non_finite_carried_state;
        return fault;
    }
};

}

#endif
