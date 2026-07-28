#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_H

/// @brief L1 adaptive controller with state predictor, projection-based adaptation,
/// and low-pass filtered control output.
///
/// @cite hovakimyan2010 -- Hovakimyan & Cao, "L1 Adaptive Control Theory", 2010, Ch. 2

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/l1_config.h"

#include "ctrlpp/dsp/vector_biquad.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <utility>

namespace ctrlpp
{

/// @brief Structured failure modes of a single `l1_controller::evaluate` cycle.
///
/// Every enumerator is an exact domain condition. The order below is the order
/// the guard tests them, and it follows the controller's own data flow, so the
/// cause named is always the most upstream one: the predictor produces the
/// prediction error, the prediction error drives the adaptation, and the
/// supplied signals enter at the end.
///
///  * non_finite_predictor_state : the state predictor's carried state is
///                                 already non-finite -- either the predicted
///                                 state or the previous filtered command, both
///                                 of which feed the predictor propagation and
///                                 are repaired the same way, by `reset`. This
///                                 is upstream of the adaptation, so it is
///                                 named ahead of it.
///  * non_finite_adaptation      : the carried uncertainty estimate is already
///                                 non-finite. Reaching this state means a NaN
///                                 survived the projection (see `evaluate` for
///                                 why an infinity would not have), so it is
///                                 specifically a NaN diagnosis and points at
///                                 the adaptation gain or the initial estimate.
///  * non_finite_state           : the supplied plant state has a non-finite
///                                 component. It enters the prediction error
///                                 and hence the adaptation. The repair is
///                                 upstream, in the sensor or the estimator.
///  * non_finite_reference       : the supplied reference command has a
///                                 non-finite component. It reaches the command
///                                 through the feedforward gain. The repair is
///                                 in whatever generates the command, a
///                                 different subsystem from the one that
///                                 measures the plant.
enum class l1_step_error
{
    non_finite_predictor_state,
    non_finite_adaptation,
    non_finite_state,
    non_finite_reference,
};

/// @brief Persistent state-health status of an `l1_controller`.
///
/// A per-cycle result cannot answer whether the controller is still carrying
/// damage, because that question outlives the call. The status latches until
/// `reset`, and it never downgrades: it is a severity ladder, so a cycle that
/// finds a lesser fault cannot mask a greater one recorded earlier.
///
///  * ok                           : every cycle so far began from finite
///                                   carried state and no adaptation update
///                                   overflowed.
///  * projection_clamped_non_finite : the adaptation update produced a
///                                   non-finite estimate and the projection
///                                   replaced it with a configured bound. This
///                                   is the quiet one and the reason this
///                                   enumerator exists at all: the cycle
///                                   SUCCEEDS, the command is finite and
///                                   in-range, and every downstream finiteness
///                                   check the caller might run passes -- while
///                                   the estimate the command was built from is
///                                   meaningless. Without this status a caller
///                                   has no way to learn it. See `evaluate`.
///  * non_finite_carried_state     : a cycle found the predictor state or the
///                                   uncertainty estimate already non-finite.
///                                   Ranked above the clamped case because the
///                                   command itself is now non-finite rather
///                                   than merely meaningless, and no cycle can
///                                   run at all until `reset`. A rejected cycle
///                                   does NOT set this by itself being
///                                   rejected on its arguments: a rejection
///                                   mutates nothing.
enum class l1_health
{
    ok,
    projection_clamped_non_finite,
    non_finite_carried_state,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX = 1, std::size_t NU = 1,
          typename Filter = vector_biquad<Scalar, NU>>
    requires vector_discrete_filter<Filter, Vector<Scalar, NU>>
class l1_controller
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NX == NU,
        "L1 reference feedforward requires NX == NU (a square predictor). The "
        "general non-square feedforward gain k_g = -(C (I - A_m)^{-1} B)^{-1} is "
        "not yet implemented");

public:
    using config_type = l1_config<Scalar, NX, NU>;
    using state_type = Vector<Scalar, NX>;
    using input_type = Vector<Scalar, NU>;

    /// Builds the output low-pass filter via `Filter::low_pass` and validates
    /// the predictor model. Returns `l1_error::invalid_filter_config` when the
    /// filter factory rejects the (cutoff_hz, sample_hz) design, or the
    /// predictor validation error of the filter-taking overload below.
    static auto create(const config_type& cfg, Scalar cutoff_hz, Scalar sample_hz)
        -> expected<l1_controller, l1_error>
    {
        auto filter = Filter::low_pass(cutoff_hz, sample_hz);
        if(!filter.has_value())
            return unexpected(l1_error::invalid_filter_config);
        return create(cfg, *std::move(filter));
    }

    /// Validates the predictor model and constructs the controller. Returns
    /// `l1_error::singular_predictor` when (I - A_m) is singular,
    /// `l1_error::singular_dc_gain` when the DC gain (I - A_m)^{-1} B is
    /// singular or non-finite, and `l1_error::non_finite_gain` when the
    /// feedforward gain K_r is non-finite.
    static auto create(const config_type& cfg, Filter filter)
        -> expected<l1_controller, l1_error>
    {
        auto k_r = compute_k_r(cfg);
        if(!k_r.has_value())
            return unexpected(k_r.error());
        return l1_controller{cfg, std::move(filter), *std::move(k_r)};
    }

    /// @brief Produce the control command for one cycle and adapt.
    ///
    /// The cycle is classified before any member is written, so a rejected
    /// cycle leaves the predictor state, the prediction error, the uncertainty
    /// estimate and the previous command bitwise unchanged, and a following
    /// valid cycle produces exactly what it would have produced had the
    /// rejected one never been attempted. As with any adaptive law, the
    /// uncertainty estimate is only accumulated into and never re-derived, so
    /// one admitted non-finite sample would be permanent.
    ///
    /// **What a refusal means for the caller.** A rejected cycle produced NO
    /// command; the actuator is still going to be driven by something and the
    /// caller must choose what. Hold the command the last successful cycle
    /// returned, command a configured safe value, or fail over -- the right
    /// answer is a property of the plant, so the controller does not pick one.
    ///
    /// **What the projection does with a non-finite estimate.** The projection
    /// is the elementwise `cwiseMax(theta_min).cwiseMin(theta_max)` below. Eigen
    /// reduces both to `numext::maxi`/`numext::mini`, which are written
    /// `x < y ? y : x` and `y < x ? y : x` (Eigen 3.4.1,
    /// `Eigen/src/Core/MathFunctions.h:1282` and `:1249`, reached through
    /// `scalar_max_op`/`scalar_min_op` at
    /// `Eigen/src/Core/functors/BinaryFunctors.h:143`). Both return the LEFT
    /// operand when the comparison is false, and every comparison against a NaN
    /// is false. Two consequences follow, and they differ:
    ///
    ///  * A NaN estimate is returned unchanged by both halves, so the
    ///    projection does NOT sanitize it. It is carried into the next cycle,
    ///    which rejects with `non_finite_adaptation`.
    ///  * An infinite estimate IS replaced -- by `theta_max` or `theta_min`,
    ///    but only when that bound is itself finite. The default configuration
    ///    leaves the bounds at -/+ infinity, and against an infinite bound the
    ///    comparison is false again and the infinity survives.
    ///
    /// The second case is the dangerous one, and it is reachable with entirely
    /// finite arguments: a large adaptation gain against a large prediction
    /// error overflows to an infinity, the projection pins it to a legitimate
    /// bound, and the controller then emits a finite, in-range command computed
    /// from an estimate that carries no information. Nothing downstream can
    /// detect that. So the raw update is tested for finiteness BEFORE the clamp
    /// and `health()` latches `projection_clamped_non_finite` when it fails,
    /// which is the only way the caller can learn it. The cycle still succeeds,
    /// because it did produce the command the algorithm prescribes; the status
    /// is what says that command should not be trusted.
    auto evaluate(const state_type& x, const input_type& r) -> expected<input_type, l1_step_error>
    {
        if(const auto step = check_step(x, r); !step)
            return unexpected(latch_health(step.error()));

        // 1. State predictor: x_hat = A_m * x_hat + B * (u_prev + sigma_hat)
        m_x_hat = propagate(m_cfg.predictor_model, m_x_hat,
                            (m_u_prev + m_sigma_hat).eval());

        // 2. Prediction error (Hovakimyan convention)
        m_x_tilde = m_x_hat - x;

        // 3. Adaptation with projection (elementwise clamp). The raw update is
        // formed first so its finiteness can be observed: the projection is
        // able to turn an infinity into a legitimate-looking bound, and that
        // substitution is invisible in every value downstream of it.
        input_type sigma_raw = m_sigma_hat;
        sigma_raw.noalias() -= m_cfg.gamma
            * (m_cfg.predictor_model.B.transpose() * m_x_tilde);
        if(!sigma_raw.allFinite())
            note_health(l1_health::projection_clamped_non_finite);
        m_sigma_hat = sigma_raw.cwiseMax(m_cfg.theta_min).cwiseMin(m_cfg.theta_max);

        // 4. Raw control: reference feedforward minus uncertainty estimate
        auto u_raw = (m_k_r * r - m_sigma_hat).eval();

        // 5. Low-pass filter (L1 robustification mechanism)
        auto u_filtered = m_filter.process(u_raw);

        // 6. Store for next predictor step
        m_u_prev = u_filtered;

        return u_filtered;
    }

    auto x_hat() const -> const state_type& { return m_x_hat; }

    auto sigma_hat() const -> const input_type& { return m_sigma_hat; }

    auto tracking_error() const -> const state_type& { return m_x_tilde; }

    auto theta() const -> const input_type& { return m_sigma_hat; }

    /// @brief Report whether the controller is still carrying damage from an
    /// earlier cycle -- either non-finite carried state, or an uncertainty
    /// estimate the projection silently substituted a bound for. Latches until
    /// `reset` and never downgrades; a cycle rejected on its arguments does not
    /// set it.
    auto health() const -> l1_health { return m_health; }

    void reset()
    {
        m_x_hat = m_cfg.x_hat_0;
        m_sigma_hat = m_cfg.sigma_hat_0;
        m_x_tilde = state_type::Zero();
        m_u_prev = input_type::Zero();
        m_filter.reset();
        // Every member the latched status describes has just been replaced, so
        // the status no longer describes anything.
        m_health = l1_health::ok;
    }

private:
    /// @brief Classify a cycle's operands without touching a single member.
    ///
    /// The order is the controller's own data flow, as documented on
    /// `l1_step_error`. The cost is one finiteness scan of each operand --
    /// 2*NX + 2*NU reads, every dimension a compile-time template parameter,
    /// with no data-dependent branching and no allocation.
    auto check_step(const state_type& x, const input_type& r) const -> expected<void, l1_step_error>
    {
        if(!m_x_hat.allFinite() || !m_u_prev.allFinite())
            return unexpected(l1_step_error::non_finite_predictor_state);
        if(!m_sigma_hat.allFinite())
            return unexpected(l1_step_error::non_finite_adaptation);
        if(!x.allFinite())
            return unexpected(l1_step_error::non_finite_state);
        if(!r.allFinite())
            return unexpected(l1_step_error::non_finite_reference);
        return {};
    }

    /// @brief Latch the persistent status for the faults that describe the
    /// controller rather than the cycle's arguments, and pass the fault through
    /// so the caller still receives the specific cause on the failure channel.
    auto latch_health(l1_step_error fault) -> l1_step_error
    {
        if(fault == l1_step_error::non_finite_predictor_state
            || fault == l1_step_error::non_finite_adaptation)
            note_health(l1_health::non_finite_carried_state);
        return fault;
    }

    /// @brief Raise the persistent status, never lower it. The enumeration is
    /// declared in severity order, so a clamped-infinity report arriving after
    /// an unrecoverable non-finite estimate cannot mask it.
    void note_health(l1_health status)
    {
        if(status > m_health)
            m_health = status;
    }

    l1_controller(const config_type& cfg, Filter filter, Matrix<Scalar, NU, NU> k_r)
        : m_cfg{cfg}
        , m_filter{std::move(filter)}
        , m_x_hat{cfg.x_hat_0}
        , m_sigma_hat{cfg.sigma_hat_0}
        , m_k_r{std::move(k_r)}
    {
    }

    static auto compute_k_r(const config_type& cfg)
        -> expected<Matrix<Scalar, NU, NU>, l1_error>
    {
        // DC gain of predictor: G_dc = (I - A_m)^{-1} * B
        // K_r = G_dc^{-1} so that in steady state x_ss = r
        auto i_minus_a = (Matrix<Scalar, NX, NX>::Identity()
            - cfg.predictor_model.A).eval();
        auto lu_ima = i_minus_a.fullPivLu();
        if(!lu_ima.isInvertible())
            return unexpected(l1_error::singular_predictor);
        Matrix<Scalar, NX, NU> dc_gain = lu_ima.solve(cfg.predictor_model.B);
        auto lu_dc = dc_gain.fullPivLu();
        if(!lu_dc.isInvertible() || !dc_gain.allFinite())
            return unexpected(l1_error::singular_dc_gain);
        Matrix<Scalar, NU, NU> k_r = lu_dc.solve(Matrix<Scalar, NU, NU>::Identity());
        if(!k_r.allFinite())
            return unexpected(l1_error::non_finite_gain);
        return k_r;
    }

    config_type m_cfg;
    Filter m_filter;
    state_type m_x_hat = state_type::Zero();
    state_type m_x_tilde = state_type::Zero();
    input_type m_sigma_hat = input_type::Zero();
    input_type m_u_prev = input_type::Zero();
    Matrix<Scalar, NU, NU> m_k_r = Matrix<Scalar, NU, NU>::Identity();
    l1_health m_health{l1_health::ok};
};

}

#endif
