#ifndef HPP_GUARD_CTRLPP_SYSID_SYSID_TYPES_H
#define HPP_GUARD_CTRLPP_SYSID_SYSID_TYPES_H

/// @brief Public error type for the batch identification routines.
///
/// `sysid_error` enumerates the structured failure modes a batch identification
/// routine can report about the data records it was handed. Together with the
/// identified model it forms the `ctrlpp::expected<arx_result<...>, sysid_error>`
/// contract of `batch_arx`.
///
/// Every enumerator is an exact precondition of the routine's index arithmetic
/// or of its data layout. None of them is a quality judgement about the fit:
/// a record that passes all four can still be uninformative, rank deficient, or
/// under-determined, and the routine will identify from it and report what it
/// found through `fit_metrics`.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

namespace ctrlpp
{

/// @brief Structured failure modes for the batch identification routines.
///
///  * record_length_mismatch : the output and input records hold different
///                             sample counts. The routine pairs sample k of one
///                             with sample k of the other, so unequal lengths
///                             have no consistent pairing.
///  * record_not_single_row  : a record has a row count other than one. These
///                             routines are single-input single-output and read
///                             only row zero, so a multi-row record would be
///                             silently identified from a fraction of its data.
///  * too_few_samples        : the sample count is at or below the larger of the
///                             two model orders. The routine forms one regressor
///                             row per sample beyond that order, so a strictly
///                             greater sample count is exactly the condition for
///                             a nonempty regressor matrix; at the order the
///                             matrix is empty and below it the row count is the
///                             result of an unsigned subtraction that wraps.
///  * non_finite_sample      : a sample in either record is NaN or infinite. The
///                             samples enter the regressor and the least-squares
///                             solve directly, which propagates the value into
///                             every identified coefficient.
enum class sysid_error
{
    record_length_mismatch,
    record_not_single_row,
    too_few_samples,
    non_finite_sample,
};

}

#endif
