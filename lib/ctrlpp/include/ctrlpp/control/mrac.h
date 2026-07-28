#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_H

/// @brief Stateful MRAC controller with Lyapunov-based adaptation and compile-time robustification.
///
/// The adaptation law projects the tracking error through B^T only, which is the
/// Lyapunov gradient B^T P e specialized to P = I. This is exact when the reference
/// model is chosen so that A_m^T + A_m is negative definite (P = I solves the
/// Lyapunov equation A_m^T P + P A_m = -Q); for a general stable A_m a different P
/// would be required, and no configurable P weighting is provided here. The update
/// also carries no explicit dt factor, so the adaptation gains gamma_x and gamma_r
/// absorb the sample time: rescale them proportionally if the sample rate changes.
///
/// @cite slotine1991 -- Slotine & Li, "Applied Nonlinear Control", 1991, Ch. 8

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/mrac_config.h"
#include "ctrlpp/control/mrac_policies.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

/// @brief Structured failure modes of a single `mrac_controller::evaluate`
/// cycle.
///
/// Every enumerator is an exact domain condition. The order below is the order
/// the guard tests them, and it follows the adaptation's own data flow, so the
/// cause named is always the most upstream one: the reference model drives the
/// tracking error, the tracking error drives both parameter matrices, and the
/// supplied signals enter at the end. Naming a downstream symptom would send
/// the caller to reset a controller whose configuration will destroy it again
/// on the next sample.
///
///  * non_finite_reference_state : the internal reference-model state is
///                                 already non-finite. It is the root cause
///                                 whenever it holds, because the tracking
///                                 error is the plant state minus this, so a
///                                 poisoned reference model poisons both
///                                 parameter matrices on the following cycle.
///                                 The repair is in the reference model or its
///                                 initial condition, not in the adaptation.
///  * non_finite_parameters      : an adaptive parameter matrix is already
///                                 non-finite. These matrices ARE the
///                                 controller's memory: nothing re-derives
///                                 them, they only accumulate, so this state is
///                                 permanent until `reset`. The repair is in
///                                 the adaptation gains or the initial
///                                 parameters.
///  * non_finite_state           : the supplied plant state has a non-finite
///                                 component. It reaches the command directly
///                                 and reaches both parameter matrices through
///                                 the tracking error, so admitting one sample
///                                 destroys the controller permanently. The
///                                 repair is upstream, in the sensor or the
///                                 estimator.
///  * non_finite_reference       : the supplied reference command has a
///                                 non-finite component. It propagates the
///                                 reference model, so it too destroys the
///                                 controller, but the repair is in whatever
///                                 generates the command -- a different
///                                 subsystem from the one that measures the
///                                 plant.
enum class mrac_step_error
{
    non_finite_reference_state,
    non_finite_parameters,
    non_finite_state,
    non_finite_reference,
};

/// @brief Persistent state-health status of an `mrac_controller`.
///
/// A per-cycle result cannot answer whether the controller is still carrying
/// damage, because that question outlives the call, and for this controller the
/// damage is the worst kind: the adaptive parameters are only ever accumulated
/// into, never re-derived, so no sequence of good samples repairs them. The
/// status latches until `reset` restores the configured initial parameters.
///
///  * ok                       : every cycle so far began from finite carried
///                               state.
///  * non_finite_carried_state : a cycle found the reference-model state or a
///                               parameter matrix already non-finite. A
///                               rejected cycle does NOT set this: a rejection
///                               mutates nothing, so it leaves the controller
///                               exactly as healthy as it was. The routes in
///                               are the non-fallible ones -- a non-finite
///                               initial condition or reference model supplied
///                               at construction -- and finite-but-extreme
///                               gains whose product overflows.
enum class mrac_health
{
    ok,
    non_finite_carried_state,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX = 1, std::size_t NU = 1,
          typename Robustification = no_robustification>
class mrac_controller
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

public:
    using config_type = mrac_config<Scalar, NX, NU, Robustification>;
    using state_type = Vector<Scalar, NX>;
    using input_type = Vector<Scalar, NU>;
    using theta_x_type = Matrix<Scalar, NU, NX>;
    using theta_r_type = Matrix<Scalar, NU, NU>;

    explicit mrac_controller(const config_type& cfg)
        : m_cfg{cfg}
        , m_x_model{cfg.x_model_0}
        , m_tracking_error{state_type::Zero()}
        , m_theta_x{cfg.theta_x_0}
        , m_theta_r{cfg.theta_r_0}
    {
    }

    /// @brief Produce the control command for one cycle and adapt.
    ///
    /// The cycle is classified before any member is written, so a rejected
    /// cycle leaves the reference-model state, the tracking error and both
    /// parameter matrices bitwise unchanged, and a following valid cycle
    /// produces exactly what it would have produced had the rejected one never
    /// been attempted.
    ///
    /// This guard matters more here than on an estimator's update. An estimator
    /// that admits one bad sample is re-driven towards the truth by the samples
    /// that follow. These parameter matrices are not: the adaptation only ever
    /// subtracts an increment into them, nothing re-derives them from data, so
    /// a single non-finite sample makes them non-finite forever and no sequence
    /// of good samples afterwards recovers. Reconstruction or `reset` is the
    /// only way back.
    ///
    /// **What a refusal means for the caller.** A rejected cycle produced NO
    /// command; the actuator is still going to be driven by something and the
    /// caller must choose what. Hold the command the last successful cycle
    /// returned, command a configured safe value, or fail over -- the right
    /// answer is a property of the plant, so the controller does not pick one.
    ///
    /// The construction path and `reset` are NOT fallible and are not converted
    /// here, so a non-finite reference model, initial parameter matrix or
    /// initial reference state still enters through them. The first cycle that
    /// follows rejects and latches `health()`, which is what that query exists
    /// for.
    auto evaluate(const state_type& x, const input_type& r) -> expected<input_type, mrac_step_error>
    {
        if(const auto step = check_step(x, r); !step)
            return unexpected(latch_health(step.error()));

        m_x_model = propagate(m_cfg.reference_model, m_x_model, r);

        m_tracking_error = x - m_x_model;

        input_type u = m_theta_x * x + m_theta_r * r;

        auto e_proj = (m_cfg.reference_model.B.transpose() * m_tracking_error).eval();

        if constexpr(std::is_same_v<Robustification, no_robustification>)
        {
            m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x;
            m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r;
        }
        else if constexpr(std::is_same_v<Robustification, dead_zone>)
        {
            if(compute_error_norm(m_tracking_error) > m_cfg.robustification.threshold)
            {
                m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x;
                m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r;
            }
        }
        else if constexpr(std::is_same_v<Robustification, sigma_modification>)
        {
            m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x
                                   + m_cfg.robustification.sigma * m_theta_x;
            m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r
                                   + m_cfg.robustification.sigma * m_theta_r;
        }
        else if constexpr(std::is_same_v<Robustification, e_modification>)
        {
            auto e_norm = compute_error_norm(m_tracking_error);
            m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x
                                   + m_cfg.robustification.delta * e_norm * m_theta_x;
            m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r
                                   + m_cfg.robustification.delta * e_norm * m_theta_r;
        }

        return u;
    }

    auto theta_x() const -> const theta_x_type& { return m_theta_x; }

    auto theta_r() const -> const theta_r_type& { return m_theta_r; }

    auto tracking_error() const -> const state_type& { return m_tracking_error; }

    auto x_model() const -> const state_type& { return m_x_model; }

    /// @brief Report whether the controller is still carrying non-finite state
    /// from an earlier cycle. Latches until `reset`; a rejected cycle does not
    /// set it.
    auto health() const -> mrac_health { return m_health; }

    void reset()
    {
        m_x_model = m_cfg.x_model_0;
        m_theta_x = m_cfg.theta_x_0;
        m_theta_r = m_cfg.theta_r_0;
        m_tracking_error = state_type::Zero();
        // Every member the latched status describes has just been replaced from
        // the configuration, so the status no longer describes anything.
        m_health = mrac_health::ok;
    }

private:
    /// @brief Classify a cycle's operands without touching a single member.
    ///
    /// The order is the adaptation's own data flow, as documented on
    /// `mrac_step_error`. The cost is one finiteness scan of each operand --
    /// NX + NU vector reads and NU*(NX + NU) matrix reads, every dimension a
    /// compile-time template parameter, with no data-dependent branching and no
    /// allocation.
    auto check_step(const state_type& x, const input_type& r) const -> expected<void, mrac_step_error>
    {
        if(!m_x_model.allFinite())
            return unexpected(mrac_step_error::non_finite_reference_state);
        if(!m_theta_x.allFinite() || !m_theta_r.allFinite())
            return unexpected(mrac_step_error::non_finite_parameters);
        if(!x.allFinite())
            return unexpected(mrac_step_error::non_finite_state);
        if(!r.allFinite())
            return unexpected(mrac_step_error::non_finite_reference);
        return {};
    }

    /// @brief Latch the persistent status for the faults that describe the
    /// controller rather than the cycle's arguments, and pass the fault through
    /// so the caller still receives the specific cause on the failure channel.
    auto latch_health(mrac_step_error fault) -> mrac_step_error
    {
        if(fault == mrac_step_error::non_finite_reference_state
            || fault == mrac_step_error::non_finite_parameters)
            m_health = mrac_health::non_finite_carried_state;
        return fault;
    }

    auto compute_error_norm(const state_type& e) const -> Scalar
    {
        return std::sqrt((e.transpose() * m_cfg.W * e)(0, 0));
    }

    config_type m_cfg;
    state_type m_x_model{};
    state_type m_tracking_error{};
    theta_x_type m_theta_x{};
    theta_r_type m_theta_r{};
    mrac_health m_health{mrac_health::ok};
};

}

#endif
