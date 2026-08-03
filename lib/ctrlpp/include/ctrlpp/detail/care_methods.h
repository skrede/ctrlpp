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
};

/// @brief Schur + Bai-Demmel reorder CARE solve path.
///
/// @note Retained for reproducibility; superseded by `sign_function_care_method`
///       after the bakeoff. The Schur path fails its primary gate by ~40
///       percent at NX=8 to 30.
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
