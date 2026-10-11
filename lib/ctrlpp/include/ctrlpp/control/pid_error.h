#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_ERROR_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_ERROR_H

/// @brief Error and status types of the PID controller.
///
/// pid_config_error is the error channel of `pid::create` and of the fallible
/// setters, pid_step_error that of a single `pid::compute` cycle, and pid_health
/// the status a controller latches across cycles.

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
///                             `reset` (or `set_integral` with a finite value
///                             for the integrator alone) recovers.
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
///                             because it is assigned into the integrator, or
///                             in the velocity form the accumulated output, so
///                             admitting it would destroy that state after an
///                             otherwise valid cycle had already committed.
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
///                              non-finite. No known public route leads here:
///                              `create`, `set_params` and `set_integral`
///                              refuse non-finite values, and a cycle whose
///                              result is non-finite commits nothing. The latch
///                              is kept as defense in depth. A rejected cycle
///                              does NOT set this: a rejection mutates nothing,
///                              so it leaves the controller exactly as healthy
///                              as it was.
enum class pid_health
{
    ok,
    non_finite_carried_state,
};

/// @brief Structured failure modes of `pid::create`, `pid::set_params` and
/// `pid::set_integral`.
///
/// The first fourteen are checked by `create` and `set_params` in the order
/// listed; a configuration with several faults reports the first. Where a
/// condition names the internal gains, those are the parallel gains: in ISA form
/// kp/Ti and kp*Td.
///
///  * non_finite_gain               : a kp, ki or kd component is NaN or
///                                    infinite. In ISA form ki is the integral
///                                    time, and a NaN or negative-infinite one
///                                    is refused; +infinity is accepted and
///                                    means no integral action (kp/Ti = 0).
///  * zero_integral_time            : in ISA form an integral time is zero. The
///                                    integral gain kp/Ti is then infinite, not
///                                    zero. No integral action is an integral
///                                    time of +infinity, or the parallel form
///                                    with ki = 0.
///  * invalid_output_limits         : an output limit is NaN, output_min
///                                    exceeds output_max in some channel,
///                                    output_min is +infinity, or output_max is
///                                    -infinity. Equal limits are accepted.
///  * non_finite_setpoint_weight    : a b or c component is NaN or infinite.
///  * invalid_filter_time_constant  : a setpoint_filter or pv_filter time
///                                    constant is NaN, negative or +infinity.
///                                    Zero means no filter.
///  * invalid_rate_limit            : a rate_limit rate_max is NaN or negative.
///                                    Zero and +infinity mean no limit.
///  * invalid_back_calculation_gain : an explicit back-calculation gain kb is
///                                    NaN, infinite or negative. Zero selects
///                                    the default gain.
///  * invalid_error_threshold       : a conditional-integration error_threshold
///                                    is not positive: NaN, zero (the default)
///                                    or negative. Such a threshold withholds
///                                    every nonzero increment, so the
///                                    integrator never moves.
///  * invalid_derivative_filter     : a deriv_filter n is NaN or infinite, or
///                                    an explicit filter time constant tf is
///                                    NaN, infinite, zero or negative.
///  * missing_oscillation_parameters : oscillation_detect is composed and an
///                                    amplitude or ultimate_frequency component
///                                    is not positive and finite. Both default
///                                    to zero, the unset value: the detector
///                                    has no scale without them.
///  * invalid_oscillation_supervision : oscillation_detect is composed and
///                                    load_count_limit or supervision_periods
///                                    is not positive and finite.
///  * undefined_tracking_time       : back_calc leaves a channel's kb at its
///                                    default, the channel's internal ki is
///                                    nonzero, and either kp and kd are both
///                                    zero or the nonzero internal kp, ki and kd
///                                    do not share a sign. The tracking time
///                                    behind the default gain then does not
///                                    exist; an explicit kb is accepted. A
///                                    channel with ki = 0 keeps kb = 0, the
///                                    limit of an infinite tracking time.
///  * undefined_derivative_filter   : deriv_filter is composed, a channel's
///                                    internal kd is nonzero and its tf is
///                                    absent, and the default Tf = kd/(kp*n) is
///                                    not a positive time: kp is zero, n is zero
///                                    or negative, or kd and kp differ in sign.
///  * unrepresentable_gain          : the configuration is well posed, and a
///                                    quantity derived from it is not
///                                    representable in Scalar: in ISA form kp
///                                    divided by the integral time, or kp times
///                                    the derivative time, overflows, or rounds
///                                    to zero while both factors are nonzero (an
///                                    integral time of +infinity excepted); the
///                                    default filter time constant kd/(kp*n) is
///                                    not a positive finite number; or the
///                                    default back-calculation gain is not
///                                    finite, or rounds to zero while the
///                                    internal ki is nonzero.
///  * non_finite_integral           : `set_integral` was handed a NaN or
///                                    infinite component.
enum class pid_config_error
{
    non_finite_gain,
    zero_integral_time,
    invalid_output_limits,
    non_finite_setpoint_weight,
    invalid_filter_time_constant,
    invalid_rate_limit,
    invalid_back_calculation_gain,
    invalid_error_threshold,
    invalid_derivative_filter,
    missing_oscillation_parameters,
    invalid_oscillation_supervision,
    undefined_tracking_time,
    undefined_derivative_filter,
    unrepresentable_gain,
    non_finite_integral,
};

}

#endif
