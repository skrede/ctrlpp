#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_METHODS_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_METHODS_H

/// @brief CARE solve-method policy tag types and concept.
///
/// Selects the algorithm used by `ctrlpp::care` and `ctrlpp::lqr_gain_continuous`
/// to solve the continuous-time algebraic Riccati equation. Three tag types
/// are provided:
///
///  * sign_function_care_method  : default, Newton iteration on sign(H) with
///                                  determinantal scaling (Roberts 1980;
///                                  Byers 1987; Higham 2008).
///  * schur_care_method          : Schur + Bai-Demmel reorder (Laub 1979;
///                                  Bai-Demmel 1993). Retained for
///                                  reproducibility; superseded by the
///                                  sign-function path after the bakeoff.
///  * balanced_schur_care_method : Schur after a DGEBAL-style diagonal
///                                  balance of the Hamiltonian (LAPACK
///                                  DGEBAL phase 2). Retained for
///                                  reproducibility; superseded by the
///                                  sign-function path after the bakeoff.
///
/// The `care_solve_method` concept constrains the trailing `Method` template
/// parameter on `care_solve_from_hamiltonian`, `care`, and `lqr_gain_continuous`
/// so overload resolution rejects foreign types.
///
/// The default assignment was picked by a governor-locked instruction-count
/// bakeoff. At the benchmarked revision, the sign-function path reached 31.8
/// to 35.8 percent fewer median instructions than `ct::optcon::CARE` at every
/// target NX in the 8 to 30 sweep; both Schur variants failed the primary gate
/// by 39 to 41 percent. The archived percentages document the method-selection
/// decision and are not a current instruction count: that revision verified
/// nothing, and the path has since gained and then partly given back work.
///
/// What the accepting path carries now, counted analytically against the same
/// cubic measure. The verification of the returned solution runs on every
/// accepted solve, where a projector-idempotence check previously let one
/// return without it. That verification is five n-by-n products at 2n^3 each,
/// three n-by-n sums at n^2, and one real Schur factorization of the closed
/// loop with the orthogonal factor NOT accumulated, about 10n^3 rather than the
/// 25n^3 accumulating it would cost (Golub and Van Loan, 4th ed., Sec. 7.5) --
/// 20n^3 in total. The projector check it replaces was one 2n-by-2n product,
/// 16n^3. The net addition is therefore 4n^3.
///
/// Against one scaled Newton iteration's own lower bound -- an LU at (2/3)m^3
/// plus a multi-right-hand-side inverse at 2m^3, with m = 2n, so (8/3)m^3 =
/// 21.33n^3 -- and the seven iterations this path was observed to take, that
/// is 4n^3 against 149.33n^3, or 2.7 percent. Accumulating the Schur factor
/// instead, as the shipped call did before it read only the quasi-triangular
/// part, would make it 19n^3 or 12.7 percent. For calibration in the same
/// accounting, the projector check that is now gone was 10.7 percent.
///
/// So the selection argument is not weakened by verifying every answer: a
/// 2.7 percent addition to the cubic work does not close a 31.8 percent margin,
/// and the alternative to paying it was a path that could report success for a
/// solution it never examined.
///
/// ## What the two Schur variants now carry, and what it did to the gap
///
/// Both Schur variants used to return their extracted matrix without examining
/// it. They now run the same verification, so each gains the full 20n^3 -- five
/// n-by-n products at 2n^3 and one non-accumulating real Schur factorization of
/// the closed loop at 10n^3 -- with nothing removed in exchange, because
/// neither had a check to replace.
///
/// Against each variant's own cubic work, counted in the same measure: the real
/// Schur factorization of the 2n-by-2n Hamiltonian WITH its orthogonal factor
/// accumulated is 25m^3 at m = 2n, or 200n^3, and the extraction's
/// rank-revealing factorization and triangular solves are 2m^3 + 3m^3 = 40n^3.
/// That is 240n^3 before the reorder is counted at all, so the addition is at
/// most 20n^3 / 240n^3 = 8.3 percent. The DGEBAL balance the second variant
/// runs first is quadratic per sweep and does not enter a cubic count.
///
/// The gap to the default therefore widened, and by how much is arithmetic
/// rather than a measurement: the default now costs 1.027 times its unverified
/// self and each Schur variant 1.083 times its own, so the ratio between them
/// grows by 1.083 / 1.027 = 1.055. The variants were already behind by the
/// archived 39 to 41 percent; they are now behind by about five and a half
/// percent more of that disadvantage. **No timing run was performed and no
/// machine-exclusivity window was requested for this** -- every number above is
/// counted, and the archived percentages remain archived rather than restated
/// as current.
///
/// ## Where a Schur variant is the better choice, measured
///
/// The instruction count is not the whole selection argument, and it used to be
/// decisive: swept over eighteen decades of a common rescale of Q and R, on a
/// comfortably damped family and on a structurally simple one, the balanced
/// variant answered every one of the 1,152 draws in each population while the
/// default answered 652 and 655 and declined the rest. **That difference no
/// longer exists.** `ctrlpp::care` now equilibrates both weightings before the
/// Hamiltonian is built, which is the same similarity the balance performs by a
/// different route, so ALL THREE TAGS ANSWER ALL 1,152 DRAWS OF EACH POPULATION.
///
/// ## What each tag equilibrates, stated by the object it rescales
///
/// The three are no longer distinguished by whether they equilibrate at all,
/// which is what the previous version of this paragraph said and what is now
/// wrong. They are distinguished by WHAT they rescale:
///
///  * `ctrlpp::care` itself, for every tag: the WEIGHTINGS, by one common
///    divisor -- their largest entry -- with the solution multiplied back
///    afterwards. That divisor is a block-diagonal similarity of the
///    Hamiltonian, so it cannot move the spectrum; see `care_weight_scale` for
///    the derivation.
///  * `balanced_schur_care_method` additionally: the HAMILTONIAN, by a
///    DGEBAL-style diagonal similarity over all 2n indices. That is strictly
///    more general than the two-block scaling above and reaches asymmetries a
///    weight divisor cannot, at the cost of an unbounded-in-data sweep whose
///    real-time properties are in the real-time matrix rather than here.
///  * the other two tags rescale nothing further.
///
/// A caller whose weights sit far from their dynamics' own scale no longer needs
/// to select a tag to get an answer. The selection argument is now the
/// instruction count above and, for the balanced variant, the sweep it runs.
///
/// Unlike the archived percentages above, THESE ACCEPTANCE COUNTS ARE PINNED
/// rather than recorded. `CARE keeps the same promise under every method tag`
/// asserts both counts as exact equalities against the drawn counts, for every
/// tag, so a tightening of the acceptance rule that costs those answers breaks a
/// test rather than leaving this paragraph wrong. It is deliberately the
/// over-rejection guard with the least margin in the suite.
///
/// Two independent criteria agree that the answers are right, and they are
/// counted separately because they establish different things. The
/// scale-invariance check compares each draw's gain against the same pose at
/// unit weight scale; it is produced by the same solver under the same tag, so
/// it cannot see an error common to both scales. The structurally simple family
/// additionally admits a closed form -- for A = -I, B = I, Q = R = sI the gain
/// is exactly (sqrt(2) - 1) I at every s -- which owes the solver nothing. Zero
/// disagreements under either criterion, for all three tags; the worst relative
/// error against the closed form is 4.020473e-16, now the same for every tag
/// because every tag solves the same equilibrated problem.
///
/// ## Reconciling the full-acceptance claim with the near-axis decline
///
/// The paragraph above says the balanced variant answers every draw of the two
/// COMMON-SCALE families. `care_error_test` says the same variant declines every
/// one of 112 poses of the NEAR-AXIS band. Both are true and they are about
/// different families. The common-scale families are well conditioned and are
/// swept by rescaling their weights, which is exactly what a diagonal balance
/// repairs. The near-axis band drives a plant eigenvalue toward the imaginary
/// axis, where the stable and unstable invariant subspaces stop being separated;
/// no balance repairs that, because the obstacle is the pose rather than its
/// scaling. A caller should read the balanced variant as answering the poses a
/// bad weight scale would otherwise lose, not as answering poses that are
/// ill-posed on their own terms.
///
/// @cite laub1979      : Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
/// @cite roberts1980   : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
/// @cite byers1987     : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
/// @cite higham2008    : Higham, "Functions of Matrices: Theory and Computation", 2008, Chapter 5
/// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 (scaling) algorithm

#include <concepts>

namespace ctrlpp::detail
{

/// @brief Default CARE solve path: Newton iteration on sign(H) with determinantal scaling.
struct sign_function_care_method
{
    /// @brief Newton steps taken before the non-contraction guard arms.
    ///
    /// After the window the iteration must not let its per-step change grow; a
    /// step that does is not converging to a sign matrix and the solve declines
    /// with `care_error::sign_function_stagnated` rather than spending the
    /// remaining budget. Inside the window the change may grow, because the
    /// determinantal scaling makes large corrections in the early steps.
    ///
    /// This is a defaulted parameter rather than a literal because the number of
    /// such early steps is a property of the input's conditioning and is not
    /// derivable from the scalar type or the dimension. Raising it trades a
    /// later decline for a chance at an answer on a badly scaled Hamiltonian;
    /// lowering it declines sooner. The guard is a bound on wasted work, not a
    /// correctness test: whatever the iteration produces is verified against the
    /// caller's own Hamiltonian before it is reported, under any value here.
    ///
    /// The default is the incumbent, and its provenance is stated rather than
    /// dressed up: it is not derived. Swept over 3,456 draws spanning eighteen
    /// decades of weight scale in each direction, every value from 0 to the
    /// iteration cap of 40 produced byte-identical outcomes -- 1,433 accepted,
    /// none unstable, 1,726 declined as stagnation and 297 through other
    /// enumerators -- so on that evidence the window has no discriminating power
    /// and 3 is retained because it moves nothing. A caller whose inputs fall
    /// outside those families is the reason this is reachable at all.
    int warmup_iterations = 3;
};

/// @brief Schur + Bai-Demmel reorder CARE solve path.
///
/// @note Retained for reproducibility; superseded by `sign_function_care_method`
///       after the bakeoff. The Schur path fails its primary gate by ~40
///       percent at NX=8 to 30. Like every other tag it verifies the matrix it
///       extracted before returning it, and reports `unverified_solution` when
///       that matrix does not satisfy the equation.
struct schur_care_method
{
};

/// @brief DGEBAL-prebalanced Schur CARE solve path.
///
/// @note Retained for reproducibility; superseded by `sign_function_care_method`
///       after the bakeoff. DGEBAL balance is a near no-op on
///       well-conditioned Hamiltonians (the diagonal D
///       scaling stays near ones, measured alongside the subspace residual), and the path
///       tracks `schur_care_method` within 1 percent across the bakeoff sweep.
///       On instruction count only, that is: on a common weight rescale far from
///       the dynamics' scale it is the only tag that answers the whole swept
///       range, as recorded above. It verifies its extracted matrix against the
///       CALLER's Hamiltonian rather than the balanced one, and reports
///       `unverified_solution` when that matrix does not satisfy the equation.
struct balanced_schur_care_method
{
};

template <typename T>
concept care_solve_method =
    std::same_as<T, schur_care_method>
 || std::same_as<T, sign_function_care_method>
 || std::same_as<T, balanced_schur_care_method>;

}

#endif
