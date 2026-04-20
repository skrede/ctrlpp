#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H

/// @brief Public types for the continuous-time algebraic Riccati equation solver.
///
/// `care_error` enumerates the structured failure modes a CARE solve can produce;
/// `care_result` carries the solution P plus diagnostic scalars. Together they form
/// the `std::expected<care_result, care_error>` contract of `ctrlpp::care`.

#include "ctrlpp/util/concepts.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp
{

/// @brief Structured failure modes for `care`.
///
///  * non_lhp_stabilisable     : fewer than n eigenvalues of the Hamiltonian spectrum
///                               lie in the open left half-plane.
///  * non_finite_input         : A, B, Q, R or the assembled Hamiltonian H contains NaN/Inf.
///  * singular_u11             : the top-left n x n block of the reordered invariant-subspace
///                               basis U is singular; P cannot be extracted.
///  * non_psd_solution         : extracted P is not positive semi-definite within an
///                               epsilon-scaled tolerance.
///  * schur_failed             : `Eigen::RealSchur` did not converge on the Hamiltonian
///                               (Schur-based methods only).
///  * sign_function_stagnated  : the matrix sign-function Newton iteration failed to
///                               contract: either the determinantal scaling factor went
///                               non-finite (singular Hamiltonian), the per-step contraction
///                               ratio exceeded 1/2 after the warm-up window (divergence),
///                               or the iteration budget was exhausted without meeting the
///                               epsilon-scaled convergence tolerance. Sign-function path only.
enum class care_error
{
    non_lhp_stabilisable,
    non_finite_input,
    singular_u11,
    non_psd_solution,
    schur_failed,
    sign_function_stagnated,
};

/// @brief Solution payload of `care`.
///
///  * P                   : n x n symmetric positive-semidefinite stabilising solution.
///  * subspace_separation : diagnostic of invariant-subspace conditioning. For Schur-based
///                          methods this is the minimum rank-revealing QR pivot ratio across
///                          all accepted block swaps during reordering (LAPACK SEP analogue);
///                          a value close to 1 indicates a well-conditioned invariant
///                          subspace, small positive values warn of near-degenerate spectra.
///                          For methods that do not run a swap phase (e.g. the matrix
///                          sign-function path), no pivot-ratio metric is defined and this
///                          field is written as `std::numeric_limits&lt;Scalar&gt;::quiet_NaN()`
///                          to signal "unavailable"; callers should branch on `std::isnan`
///                          rather than comparing against a magnitude. Contrast with the
///                          Schur path's partial-reorder marker where the smallest accepted
///                          pivot is finite and positive.
///  * reorder_complete    : true if every swap was accepted by the conditioning test or if
///                          the method has no swap phase; false if one or more swaps were
///                          declined during Schur reordering.
template <ctrlpp_floating_scalar Scalar, std::size_t NX>
struct care_result
{
    static_assert(NX > 0, "State dimension NX must be positive");

    Eigen::Matrix<Scalar, int(NX), int(NX)> P;
    Scalar                                  subspace_separation;
    bool                                    reorder_complete;
};

}

#endif
