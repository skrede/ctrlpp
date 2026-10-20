#ifndef HPP_GUARD_CTRLPP_SYSID_SYSID_TYPES_H
#define HPP_GUARD_CTRLPP_SYSID_SYSID_TYPES_H

/// @brief Public error type for the batch identification routines.
///
/// `sysid_error` enumerates the structured failure modes a batch identification
/// routine can report about the data records it was handed. Together with the
/// identified model it forms the `ctrlpp::expected<..., sysid_error>` contract
/// of `batch_arx`, `moesp` and `moesp_singular_values`.
///
/// Every enumerator except `model_not_identified` is an exact precondition of
/// the routine's index arithmetic or of its data layout. None of those is a
/// quality judgement about the fit: a record that passes them can still be
/// uninformative, rank deficient, or under-determined, and the routine will
/// identify from it and report what it found through `fit_metrics`.
/// `model_not_identified` reports that the identification itself produced no
/// model of the requested order.
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
///  * too_few_samples        : for `batch_arx`, the sample count is at or below
///                             the larger of the two model orders. The routine
///                             forms one regressor row per sample beyond that
///                             order, so a strictly greater sample count is
///                             exactly the condition for a nonempty regressor
///                             matrix; at the order the matrix is empty and below
///                             it the row count is the result of an unsigned
///                             subtraction that wraps. For `moesp`, the record
///                             cannot fill the block Hankel matrices the
///                             block-row count i asks for with the NX projection
///                             columns the order needs: N - 2i + 1 < i + NX,
///                             with NX = 1 for `moesp_singular_values`. At
///                             N - 2i + 1 <= 0 the matrices do not exist at all.
///  * non_finite_sample      : a sample in either record is NaN or infinite. The
///                             samples enter the regressor and the least-squares
///                             solve directly, which propagates the value into
///                             every identified coefficient.
///  * block_rows_below_order : the shift relation Gamma_trunc A = Gamma_shifted
///                             has (i - 1) * ny rows and NX unknowns per column,
///                             so below NX rows it is underdetermined and its
///                             minimum-norm solution is finite and meaningless.
///                             The default block-row count min(N / 4, 30) is
///                             checked too, so a short record without an
///                             explicit count can report this.
///  * model_not_identified   : an extracted system matrix is non-finite, or the
///                             condition number sv(0) / sv(NX - 1) of the
///                             projection is not finite: its NX-th singular value
///                             is zero, or so small against the leading one that
///                             the ratio overflows, and the NX-th direction of
///                             the subspace is undetermined.
enum class sysid_error
{
    record_length_mismatch,
    record_not_single_row,
    too_few_samples,
    non_finite_sample,
    block_rows_below_order,
    model_not_identified,
};

}

#endif
