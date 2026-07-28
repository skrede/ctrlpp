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
/// bakeoff. The sign-function path reached 31.8 to 35.8 percent fewer median
/// instructions than `ct::optcon::CARE` at every target NX in the 8 to 30
/// sweep; both
/// Schur variants failed the primary gate by 39 to 41 percent.
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
