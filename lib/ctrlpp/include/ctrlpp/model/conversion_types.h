#ifndef HPP_GUARD_CTRLPP_MODEL_CONVERSION_TYPES_H
#define HPP_GUARD_CTRLPP_MODEL_CONVERSION_TYPES_H

/// @brief Public types for the transfer-function and state-space conversions.
///
/// `conversion_error` enumerates the structured refusals `ctrlpp::tf2ss` and
/// `ctrlpp::ss2tf` can produce, one enumerator per distinct cause. Together with
/// the routines it forms their `ctrlpp::expected<..., conversion_error>` contract.

namespace ctrlpp
{

/// @brief Structured refusals for `tf2ss` and `ss2tf`.
///
///  * non_finite_input          : a transfer-function coefficient or a
///                                state-space entry is NaN or infinite.
///  * zero_leading_coefficient  : the denominator's leading coefficient is
///                                zero, so its true degree is below the
///                                realization's dimension and the companion
///                                form of that dimension does not exist.
///  * eigenvalues_not_converged : the eigensolver did not converge, so no
///                                characteristic polynomial was formed.
///  * unrepresentable_result    : the input is finite but a converted
///                                coefficient is not (a tiny leading
///                                coefficient's quotient, or a polynomial
///                                coefficient beyond the scalar's range).
enum class conversion_error
{
    non_finite_input,
    zero_leading_coefficient,
    eigenvalues_not_converged,
    unrepresentable_result,
};

}

#endif
