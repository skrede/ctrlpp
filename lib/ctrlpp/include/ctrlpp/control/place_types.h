#ifndef HPP_GUARD_CTRLPP_CONTROL_PLACE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_PLACE_TYPES_H

/// @brief Public types for the pole-placement routines.
///
/// `place_error` enumerates the structured refusals `ctrlpp::place` and
/// `ctrlpp::place_observer` can produce, one enumerator per distinct cause.
/// Together with the routines it forms their
/// `ctrlpp::expected<Eigen::Matrix<...>, place_error>` contract.

namespace ctrlpp
{

/// @brief Structured refusals for `place` and `place_observer`.
///
/// The first two are STRUCTURAL: the caller cannot reach a gain by changing the
/// numbers, because the single-channel Ackermann formula the routines implement
/// does not cover the shape at all. The last two are conditions ON THE DATA,
/// which a different pole set or a different pair can satisfy. Sharing an
/// enumerator across that boundary would send the caller to change the one
/// thing that cannot help.
///
///  * multi_input_not_supported     : `place` was instantiated with an input
///                                    dimension above one. Ackermann's formula
///                                    assigns a characteristic polynomial
///                                    through a single input channel, and the
///                                    multi-input assignment is not a unique
///                                    problem -- it has a free subspace, which
///                                    is what a robust-assignment method
///                                    exists to choose within.
///  * multi_output_not_supported    : `place_observer` was instantiated with an
///                                    output dimension above one. The dual of
///                                    the above, through the same formula.
///  * poles_not_conjugate_symmetric : the requested pole set is not closed
///                                    under conjugation, so no real-coefficient
///                                    characteristic polynomial has it as its
///                                    root set and no real gain assigns it.
///  * uncontrollable_pair           : the controllability matrix
///                                    [B, AB, ..., A^{n-1} B] is rank-deficient
///                                    to its rank-revealing QR threshold, so
///                                    the transform Ackermann's formula inverts
///                                    does not exist and some modes cannot be
///                                    moved at all. For `place_observer` the
///                                    same enumerator reports the dual
///                                    condition on (A^T, C^T), which is
///                                    unobservability.
enum class place_error
{
    multi_input_not_supported,
    multi_output_not_supported,
    poles_not_conjugate_symmetric,
    uncontrollable_pair,
};

}

#endif
