#ifndef HPP_GUARD_CTRLPP_CONTROL_PID_ERROR_H
#define HPP_GUARD_CTRLPP_CONTROL_PID_ERROR_H

namespace ctrlpp
{

// A cycle tests the first four in declaration order, the most upstream cause first: the carried state,
// then the clock, then the two signals.
enum class pid_step_error
{
    non_finite_state,
    invalid_timestep,
    non_finite_setpoint,
    non_finite_measurement,
    non_finite_tracking_signal,
    non_finite_result,
};

// Latches until reset. A rejected cycle mutates nothing and so does not set it.
enum class pid_health
{
    ok,
    non_finite_carried_state,
};

// create and set_params test the first fourteen in declaration order and report the first that fails.
// A condition on a gain reads the internal parallel gains: in ISA form kp/Ti and kp*Td.
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
