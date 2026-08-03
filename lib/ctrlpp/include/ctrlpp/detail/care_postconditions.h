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

#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cstddef>

namespace ctrlpp::detail
{

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
    // residual evaluation. A conservative count uses 2n^3 operations for the
    // n-by-n Householder QR, 3n^3 for applying Q and solving n triangular
    // right-hand sides, and 3n^2 for symmetrization. Four n-by-n products then
    // use n^2 * (2n - 1) operations each, and three n-by-n sums combine the
    // residual terms. Every operation is counted whether or not it rounds.
    // The Frobenius scale is the largest term entering the final sum, not the
    // cancellation residual.
    constexpr int extraction_qr_rounding_ops =
        2 * n * n * n;
    constexpr int extraction_solve_rounding_ops =
        3 * n * n * n;
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
        extraction_qr_rounding_ops
        + extraction_solve_rounding_ops
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
    // quasi-triangular: a nonzero subdiagonal entry marks a 2 x 2 block holding
    // a complex conjugate pair, and Eigen does not standardize such a block, so
    // its two diagonal entries are not the pair's real part -- only their mean
    // is. Reading the diagonal directly refuses every oscillatory closed loop,
    // including the exact solution of the double integrator, whose factor
    // carries the diagonal (0, -sqrt(3)) for a pair whose real parts are both
    // -sqrt(3)/2. The blocks are therefore walked rather than the diagonal.
    constexpr int closed_loop_eigenvalue_rounding_ops = n;
    const resolved_magnitude<Scalar> closed_loop_margin =
        counted_rounding_bound(resolve_largest_entry(closed_loop_T),
                               closed_loop_eigenvalue_rounding_ops);
    if (!closed_loop_margin.resolved)
        return false;

    int index = 0;
    while (index < n)
    {
        const bool is_conjugate_pair =
            index + 1 < n
            && closed_loop_T(index + 1, index) != Scalar{0};
        const Scalar eigenvalue_real_part =
            is_conjugate_pair
                ? (closed_loop_T(index, index)
                   + closed_loop_T(index + 1, index + 1))
                      / Scalar{2}
                : closed_loop_T(index, index);
        if (!(eigenvalue_real_part < -closed_loop_margin.value))
            return false;
        index += is_conjugate_pair ? 2 : 1;
    }
    return true;
}

}

#endif
