#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_TYPES_H

/// @brief Public types for the discrete algebraic Riccati equation solver.
///
/// `dare_error` enumerates the structured failure modes a DARE solve can produce;
/// `dare_result` carries the solution P plus diagnostic scalars (subspace separation,
/// reorder completeness). Together they form the `std::expected<dare_result, dare_error>`
/// contract of `ctrlpp::dare`.

#include <Eigen/Core>

#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

/// @brief Structured failure modes for `dare`.
///
///  * non_stabilisable : fewer than n eigenvalues of the symplectic spectrum lie in
///                       the stable (|lambda| < 1) region.
///  * non_finite_input : A, B, Q, R or the assembled symplectic Z contains NaN/Inf,
///                       or A is singular at symplectic build time.
///  * singular_u11     : the top-left n x n block of the reordered invariant-subspace
///                       basis U is singular; P cannot be extracted.
///  * non_psd_solution : extracted P is not positive semi-definite within an
///                       epsilon-scaled tolerance.
///  * schur_failed     : `Eigen::RealSchur` did not converge on the symplectic matrix.
enum class dare_error
{
    non_stabilisable,
    non_finite_input,
    singular_u11,
    non_psd_solution,
    schur_failed,
};

/// @brief Solution payload of `dare`. Carries the Riccati solution P plus diagnostic
/// scalars mirroring LAPACK DTRSEN's `SEP` and `INFO=1` semantics.
///
///  * P                   : n x n symmetric positive-semidefinite stabilising solution.
///  * subspace_separation : minimum rank-revealing QR pivot ratio across all accepted
///                          block swaps during Schur reordering (LAPACK SEP analogue).
///                          A value close to 1 indicates a well-conditioned invariant
///                          subspace; small values warn of near-degenerate spectra.
///  * reorder_complete    : true if every swap was accepted by the conditioning test;
///                          false if one or more swaps were declined (LAPACK INFO=1
///                          analogue). A partial reorder with a computable P is
///                          diagnostic, not an error; consult `subspace_separation`
///                          to decide whether P is trustworthy for the use case.
template <typename Scalar, std::size_t NX>
struct dare_result
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");

    Eigen::Matrix<Scalar, int(NX), int(NX)> P;
    Scalar                                  subspace_separation;
    bool                                    reorder_complete;
};

}

#endif
