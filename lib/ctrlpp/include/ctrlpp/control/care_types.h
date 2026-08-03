#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H

/// @brief Public types for the continuous-time algebraic Riccati equation solver.
///
/// `care_error` enumerates the structured failure modes a CARE solve can produce;
/// `care_result` carries the solution P plus diagnostic scalars. Together they form
/// the `ctrlpp::expected<care_result, care_error>` contract of `ctrlpp::care`.

#include "ctrlpp/util/concepts.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp
{

/// @brief Structured failure modes for `care`.
///
///  * non_lhp_stabilizable     : fewer than n eigenvalues of the Hamiltonian spectrum
///                               lie in the open left half-plane. See the note below on
///                               what this test can and cannot see.
///  * non_finite_input         : A, B, Q, R or the assembled Hamiltonian H contains NaN/Inf.
///  * singular_r               : the input weighting R is rank-deficient to a
///                               scale-relative reciprocal-pivot tolerance, so the
///                               R^{-1} the Hamiltonian build requires for
///                               B R^{-1} B^T does not exist. A caller who set a
///                               weighting to zero deliberately should not be told
///                               their input was non-finite, nor handed a solution to
///                               a different problem.
///  * singular_u11             : the top-left n x n block of the reordered invariant-subspace
///                               basis U is singular; P cannot be extracted. See the
///                               note below: this covers two different situations and
///                               does not distinguish them.
///  * non_psd_solution         : extracted P is not positive semi-definite within an
///                               epsilon-scaled tolerance.
///  * schur_failed             : `Eigen::RealSchur` did not converge on the Hamiltonian
///                               (Schur-based methods only).
///  * sign_function_stagnated  : the matrix sign-function Newton iteration did not
///                               produce a solution the solver could verify. Exactly
///                               four cases reach it: the determinantal scaling factor
///                               went non-finite (singular Hamiltonian); the per-step
///                               contraction ratio exceeded 1/2 after the warm-up window
///                               (divergence); an acceptance magnitude could not be
///                               resolved at the scalar type's range, so the comparison
///                               that would have used it carries no evidence; or the
///                               extracted solution failed the counted Riccati-residual
///                               or closed-loop-spectrum postconditions. Those
///                               postconditions run on every solve, so the last case now
///                               covers every unverified answer rather than only those a
///                               projector check had left unresolved. The iteration
///                               budget being exhausted also reports here.
///                               **Sign-function path only, and deliberately so**: every
///                               one of its four cases is a statement about a Newton
///                               iteration, so a Schur extraction that failed the same
///                               postconditions reports `unverified_solution` instead
///                               rather than being described as an iteration that
///                               stagnated.
///  * unverified_solution      : the method extracted a candidate solution, and that
///                               matrix failed the postconditions every continuous method
///                               is held to -- the counted Riccati residual bound, the
///                               requirement that the closed-loop spectrum lie strictly in
///                               the open left half-plane, or the resolvability of an
///                               acceptance magnitude at the scalar type's range.
///                               Produced by `schur_care_method` and
///                               `balanced_schur_care_method`. It says nothing about which
///                               algorithm ran and is NOT the sign-function stagnation
///                               cause: it reports that an answer was produced and could
///                               not be verified, not that an iteration failed to
///                               converge. A caller who meets it has selected a method
///                               whose extraction did not resolve this problem accurately
///                               enough to certify; the default sign-function tag and the
///                               balanced variant are the alternatives to try.
///
/// ## What `singular_u11` covers, and what `non_lhp_stabilizable` misses
///
/// The continuous solver has the identical shape as its discrete counterpart, for the
/// identical reason. `non_lhp_stabilizable` fires when fewer than n eigenvalues of the
/// Hamiltonian spectrum lie in the open left half-plane. An uncontrollable mode at
/// Re(lambda) > 0 contributes BOTH lambda and its reflection -lambda, and the
/// reflection IS in the left half-plane, so the count is satisfied and this enumerator
/// never fires for that pair. What fails instead is the extraction, and the refusal
/// arrives as `singular_u11`. Measured: A = diag(2, -1/2), B = [0; 1], Q = I, R = 1 is
/// uncontrollable in its unstable mode and refuses with `singular_u11`.
///
/// The pair is still refused, so this is a naming limit rather than a correctness one.
/// The implication that makes `singular_u11` informative is exact: a stabilizable and
/// detectable pair has a nonsingular U11 (Laub 1979 Sec. III), so **in exact
/// arithmetic** a singular U11 implies the pair is not both stabilizable and
/// detectable.
///
/// **That qualifier is load-bearing, and the enumerator is therefore NOT renamed.**
/// The test is a numerical rank test with a threshold relative to the largest pivot,
/// so a genuinely well-posed pair whose invariant subspace is severely ill-conditioned
/// reaches the same branch. `singular_u11` covers both situations and **does not
/// distinguish them**; telling them apart needs a stabilizability test the solver does
/// not perform.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979, Sec. III
enum class care_error
{
    non_lhp_stabilizable,
    non_finite_input,
    singular_r,
    singular_u11,
    non_psd_solution,
    schur_failed,
    sign_function_stagnated,
    unverified_solution,
};

/// @brief Solution payload of `care`.
///
/// ## The promise `P` carries, and which methods keep it
///
/// `P` is the stabilizing solution, and that statement is **not qualified by the
/// method tag**. Every selectable method verifies the matrix it extracted, against
/// the caller's own Hamiltonian, before reporting success: the same counted Riccati
/// residual bound and the same strict open-left-half-plane requirement on the closed
/// loop, from one definition in `ctrlpp/detail/care_postconditions.h`. A method that
/// cannot verify what it extracted declines instead of returning it.
///
/// This was previously true of the default sign-function tag alone. The two Schur
/// tags returned their extracted matrix unverified, and on a controllable and
/// detectable near-axis family they returned an unstable closed loop in 150 of 814
/// and 156 of 839 successes. A caller who selected a method tag was therefore
/// selecting a different contract without being told so; they no longer are.
///
///  * P                   : n x n symmetric positive-semidefinite stabilizing solution.
///  * subspace_separation : diagnostic of invariant-subspace conditioning. For Schur-based
///                          methods this is the minimum rank-revealing QR pivot ratio across
///                          all accepted block swaps during reordering (LAPACK SEP analogue);
///                          a value close to 1 indicates a well-conditioned invariant
///                          subspace, small positive values warn of near-degenerate spectra.
///                          For methods that do not run a swap phase (e.g. the matrix
///                          sign-function path), no pivot-ratio metric is defined and this
///                          field is written as `std::numeric_limits&lt;Scalar&gt;::quiet_NaN()`
///                          to signal "unavailable"; callers should branch on `std::isnan`
///                          rather than comparing against a magnitude. The sign path's
///                          acceptance decision is a verification of the returned solution
///                          rather than a conditioning measurement, so it yields no value
///                          for this field either. Contrast with the Schur path's
///                          partial-reorder marker where the smallest accepted pivot is
///                          finite and positive.
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
