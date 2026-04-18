#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_METHODS_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_METHODS_H

/// @brief CARE solve-method policy tag types and concept.
///
/// Selects the algorithm used by `ctrlpp::care` and `ctrlpp::lqr_gain_continuous`
/// to solve the continuous-time algebraic Riccati equation. Three tag types
/// are provided:
///
///  * schur_care_method          : default, Schur + Bai-Demmel reorder
///                                  (Laub 1979; Bai-Demmel 1993); the
///                                  post-integration baseline.
///  * sign_function_care_method  : Newton iteration on sign(H) with
///                                  determinantal scaling (Roberts 1980;
///                                  Byers 1987; Higham 2008).
///  * balanced_schur_care_method : Schur after a DGEBAL-style diagonal
///                                  balance of the Hamiltonian (LAPACK
///                                  DGEBAL phase 2).
///
/// The `care_solve_method` concept constrains the trailing `Method` template
/// parameter on `care_solve_from_hamiltonian`, `care`, and `lqr_gain_continuous`
/// so overload resolution rejects foreign types.
///
/// @cite laub1979      : Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
/// @cite roberts1980   : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
/// @cite byers1987     : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
/// @cite higham2008    : Higham, "Functions of Matrices: Theory and Computation", 2008, Chapter 5
/// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 (scaling) algorithm

#include <concepts>
#include <type_traits>

namespace ctrlpp::detail
{

struct schur_care_method
{
};

struct sign_function_care_method
{
};

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
