#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_TYPES_H

/// @brief Public types for the discrete algebraic Riccati equation solver.
///
/// `dare_error` enumerates the structured failure modes a DARE solve can produce;
/// `dare_result` carries the solution P plus diagnostic scalars (subspace separation,
/// reorder completeness). Together they form the `ctrlpp::expected<dare_result, dare_error>`
/// contract of `ctrlpp::dare`.

#include "ctrlpp/util/concepts.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp {

/// @brief Structured failure modes for `dare`.
///
///  * non_stabilisable : fewer than n eigenvalues of the symplectic spectrum lie in
///                       the stable (|lambda| < 1) region. See the note below on what
///                       this test can and cannot see.
///  * non_finite_input : A, B, Q or R contains NaN/Inf.
///  * singular_a       : the state matrix A is rank-deficient to a scale-relative
///                       reciprocal-pivot tolerance, so the A^{-T} the symplectic
///                       pencil build requires (Laub Eq. 7) does not exist.
///  * singular_r       : the input weighting R is rank-deficient to the same
///                       scale-relative reciprocal-pivot tolerance, so the R^{-1} the
///                       same pencil build requires for G = B R^{-1} B^T does not
///                       exist. The exact analogue of singular_a on the other
///                       inverted operand, and it exists because a caller who set a
///                       weighting to zero deliberately should not be told their
///                       input was non-finite.
///  * singular_u11     : the top-left n x n block of the reordered invariant-subspace
///                       basis U is singular; P cannot be extracted. See the note
///                       below: this covers two different situations and does not
///                       distinguish them.
///  * non_psd_solution : extracted P is not positive semi-definite within an
///                       epsilon-scaled tolerance.
///  * schur_failed     : `Eigen::RealSchur` did not converge on the symplectic matrix.
///  * arithmetic_limit : a finite, structurally admissible problem could not produce
///                       a solution whose Riccati residual and stabilizing closed-loop
///                       spectrum are defensible at the scalar type's precision.
///
/// ## What `singular_u11` covers, and what `non_stabilisable` misses
///
/// These two are worth reading together, because the placement count cannot see
/// every pair that has no stabilising solution.
///
/// `non_stabilisable` fires when fewer than n eigenvalues of the symplectic spectrum
/// lie inside the unit disk. An uncontrollable mode at |lambda| = 1 contributes two
/// eigenvalues ON the circle, neither inside, so the count falls short and the
/// enumerator fires. But an uncontrollable mode at |lambda| > 1 contributes BOTH
/// lambda and its reciprocal, and the reciprocal IS inside the disk, so the count is
/// satisfied and this enumerator never fires for that pair. What fails instead is the
/// extraction: the invariant subspace does not project onto the state space, and the
/// refusal arrives as `singular_u11`.
///
/// The pair is still refused -- no gain is ever returned for a pair with no
/// stabilising solution -- so this is a naming limit, not a correctness one. The
/// implication that makes `singular_u11` informative is standard and exact: a
/// stabilisable and detectable pair has a nonsingular U11 (Laub 1979 Sec. III), so
/// **in exact arithmetic** a singular U11 implies the pair is not both stabilisable
/// and detectable.
///
/// **That qualifier is load-bearing, and the enumerator is therefore NOT renamed.**
/// The test performed here is a numerical rank test with a threshold relative to the
/// largest pivot, so a genuinely well-posed pair whose invariant subspace is severely
/// ill-conditioned reaches the same branch. Measured on the one-parameter family
/// A = diag(2, 1/2), B = [delta; 1], Q = I, R = 1, which is controllable and
/// detectable for every delta != 0: the solve succeeds down to delta = 1e-7 with
/// ||P|| ~ 8e14, and returns `singular_u11` from delta = 1e-8 downward, where the
/// smallest pivot of the computed U11 falls under the rank test's threshold. Naming
/// this enumerator after the structural cause would state something false about every
/// one of those rows.
///
/// So: `singular_u11` covers a pair with no stabilising solution AND a numerical
/// failure to separate the invariant subspace, and **it does not distinguish them**.
/// A caller that must tell them apart needs a stabilisability test the solver does not
/// perform.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979, Sec. III
enum class dare_error
{
    non_stabilisable,
    non_finite_input,
    singular_a,
    singular_r,
    singular_u11,
    non_psd_solution,
    schur_failed,
    arithmetic_limit,
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
template<ctrlpp_floating_scalar Scalar, std::size_t NX>
struct dare_result
{
    static_assert(NX > 0, "State dimension NX must be positive");

    Eigen::Matrix<Scalar, int(NX), int(NX)> P;
    Scalar subspace_separation;
    bool reorder_complete;
};

}

#endif
