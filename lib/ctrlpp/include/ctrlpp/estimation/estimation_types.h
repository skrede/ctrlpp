#ifndef HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H

namespace ctrlpp
{

enum class filter_error
{
    degenerate_quaternion,
    non_positive_sigma_spread,
    non_positive_scaling_radicand,
    non_finite_process_noise,
    non_finite_measurement_noise,
    non_finite_initial_state,
    non_finite_initial_covariance,
    indefinite_process_noise,
    indefinite_measurement_noise,
    singular_measurement_noise,
    indefinite_initial_covariance,
    asymmetric_process_noise,
    asymmetric_measurement_noise,
    asymmetric_initial_covariance,
    non_unit_quaternion,
    non_finite_initial_bias,
};

}

#endif
