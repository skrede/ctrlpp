#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_POSTCONDITIONS_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_POSTCONDITIONS_H

/// @brief THE acceptance rule for a continuous Riccati solution, defined once.
///
/// A solution is accepted when it satisfies the equation it was asked to solve
/// and places the closed loop strictly in the open left half-plane. Both are
/// properties of the RETURNED MATRIX, evaluated against the caller's own
/// Hamiltonian, so the rule says nothing about which algorithm produced the
/// matrix and cannot be satisfied by an inference about that algorithm's
/// internal state.
///
/// That method-independence is why the rule lives in its own header. Every
/// continuous solve path reaches it, and none of them reaches it by including
/// another path's implementation: the sign-function path, the real-Schur path
/// and the balanced-Schur path each include this file directly. A second copy
/// specialized to one path is the failure this arrangement exists to prevent --
/// the published result contract promises a stabilizing matrix without
/// qualifying by method, and two copies are two contracts.
///
/// Both legs are counted rather than tuned. The residual is compared against
/// an enumerated operation count times unit roundoff times the largest term
/// entering the residual sum, and the closed-loop margin against the same
/// counted form at the Schur factor's own scale. Every magnitude is a
/// ctrlpp::detail::resolved_magnitude, so neither end of the arithmetic range
/// can turn a comparison into an unconditional verdict: an unresolved magnitude
/// declines rather than certifies.
///
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.5

#include "ctrlpp/detail/quasi_triangular.h"
#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cstddef>

namespace ctrlpp::detail
{

/// @brief Dimension of the factorization whose orthogonal factor becomes the
/// basis every continuous extraction reads.
///
/// It is 2n on every path, not n. The real-Schur and balanced-Schur paths
/// factorize the 2n-by-2n Hamiltonian and hand over its leading n columns; the
/// sign-function path factorizes the 2n-by-2n projector (I - sign(H))/2, whose
/// orthogonal factor IS that basis. Naming it once is the point: a bound
/// written at n-by-n under-counts every path by roughly a factor of eight, and
/// the same count is quoted by the convergence anchor rather than respelled
/// there with its own basis size.
template <std::size_t NX>
inline constexpr int care_extraction_basis_dimension = 2 * int(NX);

/// @brief Rounded operations in the extraction that produced the continuous
/// Riccati solution: 2m^3 for the m-by-m Householder QR and 3m^3 for applying
/// the orthogonal factor and solving the triangular right-hand sides, at
/// m = care_extraction_basis_dimension.
template <std::size_t NX>
inline constexpr int care_extraction_rounding_ops =
    2 * care_extraction_basis_dimension<NX>
      * care_extraction_basis_dimension<NX>
      * care_extraction_basis_dimension<NX>
    + 3 * care_extraction_basis_dimension<NX>
        * care_extraction_basis_dimension<NX>
        * care_extraction_basis_dimension<NX>;

template <typename Scalar, std::size_t NX>
auto care_solution_satisfies_postconditions(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H,
    const Eigen::Matrix<Scalar, int(NX), int(NX)>& P) -> bool
{
    constexpr int n = int(NX);
    using MatN = Eigen::Matrix<Scalar, n, n>;

    const MatN A = H.template block<n, n>(0, 0);
    const MatN S = -H.template block<n, n>(0, n);
    const MatN Q = -H.template block<n, n>(n, 0);
    const MatN At_P = (A.transpose() * P).eval();
    const MatN P_A = (P * A).eval();
    const MatN P_S_P = (P * S * P).eval();
    const MatN residual = (At_P + P_A - P_S_P + Q).eval();

    // The extracted P carries the arithmetic of its solve as well as the
    // residual evaluation. The extraction count is the shared one above, taken
    // at the 2n basis every path's factorization actually produces. On top of
    // it: 3n^2 for symmetrization, four n-by-n products at n^2 * (2n - 1)
    // operations each, and three n-by-n sums combining the residual terms.
    // Every operation is counted whether or not it rounds. The Frobenius scale
    // is the largest term entering the final sum, not the cancellation
    // residual.
    constexpr int extraction_rounding_ops = care_extraction_rounding_ops<NX>;
    constexpr int symmetrize_rounding_ops = 3 * n * n;
    constexpr int matrix_product_entry_rounding_ops = 2 * n - 1;
    constexpr int residual_matrix_products = 4;
    constexpr int residual_product_rounding_ops =
        residual_matrix_products * n * n
        * matrix_product_entry_rounding_ops;
    constexpr int residual_sums = 3;
    constexpr int residual_sum_rounding_ops =
        residual_sums * n * n;
    constexpr int residual_rounding_ops =
        extraction_rounding_ops
        + symmetrize_rounding_ops
        + residual_product_rounding_ops
        + residual_sum_rounding_ops;
    const resolved_magnitude<Scalar> residual_scale =
        largest_magnitude<Scalar>({resolve_magnitude(At_P),
                                   resolve_magnitude(P_A),
                                   resolve_magnitude(P_S_P),
                                   resolve_magnitude(Q)});
    const resolved_magnitude<Scalar> residual_floor =
        counted_rounding_bound(residual_scale, residual_rounding_ops);
    if (!magnitude_within(resolve_magnitude(residual), residual_floor))
        return false;

    // Only the quasi-triangular factor is read below, so the orthogonal factor
    // is not accumulated. Eigen's constructor defaults to computing it, which
    // costs about 25n^3 against the 10n^3 of the factor alone.
    const MatN closed_loop = (A - S * P).eval();
    Eigen::RealSchur<MatN> closed_loop_schur(closed_loop, false);
    if (closed_loop_schur.info() != Eigen::Success)
        return false;

    const MatN closed_loop_T = closed_loop_schur.matrixT();
    if (!closed_loop_T.allFinite())
        return false;

    // A backward-stable real Schur factor has an n-operation eigenvalue
    // uncertainty at the factor's largest-entry scale. Requiring every
    // eigenvalue's real part below its negative bound certifies that the
    // closed-loop spectrum lies strictly in the open left half-plane.
    //
    // The diagonal of a real Schur factor is NOT that spectrum. The factor is
    // quasi-triangular, and a 2 x 2 block's eigenvalues are its half-trace plus
    // and minus the square root of its discriminant. Reading the diagonal
    // directly refuses every oscillatory closed loop, including the exact
    // solution of the double integrator, whose factor carries the diagonal
    // (0, -sqrt(3)) for a pair whose real parts are both -sqrt(3)/2. Reading
    // the half-trace alone is wrong in the opposite and more dangerous
    // direction: a block whose discriminant is non-negative holds two REAL
    // roots, and one of them can sit in the open right half-plane while their
    // mean sits left of any margin. The walk therefore asks the discriminant
    // and tests the root with the larger real part, from the shared primitive
    // in detail/quasi_triangular.h that the reordering path also uses.
    constexpr int closed_loop_eigenvalue_rounding_ops = n;
    const resolved_magnitude<Scalar> closed_loop_margin =
        counted_rounding_bound(resolve_largest_entry(closed_loop_T),
                               closed_loop_eigenvalue_rounding_ops);
    if (!closed_loop_margin.resolved)
        return false;

    return quasi_triangular_spectrum_strictly_left_of<Scalar, n>(
        closed_loop_T, closed_loop_margin.value);
}

}

#endif
