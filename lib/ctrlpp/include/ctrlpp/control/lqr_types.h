#ifndef HPP_GUARD_CTRLPP_CONTROL_LQR_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_LQR_TYPES_H

/// @brief Public types for the linear-quadratic sequence helpers.
///
/// `lqr_error` enumerates the structured refusals `ctrlpp::lqr_tv_gains`,
/// `ctrlpp::lqr_cost`, `ctrlpp::lqr_time_varying::compute` and
/// `ctrlpp::lqr_time_varying::gain` can produce. Together with those helpers it
/// forms their `ctrlpp::expected<..., lqr_error>` contract. Every enumerator is
/// an exact precondition of the helpers' indexing or of their arithmetic's
/// domain, so none carries a tolerance.

namespace ctrlpp
{

/// @brief Structured refusals for the linear-quadratic sequence helpers.
///
///  * sequence_length_mismatch   : the four per-step sequences passed to
///                                 `lqr_tv_gains` hold different counts, so
///                                 some step has no consistent set of system
///                                 and weight matrices.
///  * trajectory_length_mismatch : the state sequence passed to `lqr_cost`
///                                 holds neither as many entries as the input
///                                 sequence nor exactly one more.
///  * step_out_of_range          : the step passed to `lqr_time_varying` is at
///                                 or past the stored horizon.
///  * non_finite_input           : a matrix or vector entry is NaN or
///                                 infinite.
enum class lqr_error
{
    sequence_length_mismatch,
    trajectory_length_mismatch,
    step_out_of_range,
    non_finite_input,
};

}

#endif
