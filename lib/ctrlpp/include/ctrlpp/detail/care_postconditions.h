#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_POSTCONDITIONS_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_POSTCONDITIONS_H

// The acceptance rule every continuous Riccati path applies to its returned matrix, against the caller's own
// Hamiltonian: the residual within a counted rounding bound, and the closed-loop spectrum strictly in the open left
// half-plane. It is kept apart from every path so the result contract has one definition whatever the method.
// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.5

#include "ctrlpp/detail/quasi_triangular.h"
#include "ctrlpp/detail/riccati_solution.h"

#include "ctrlpp/control/care_types.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

// 2n on every path: the Schur paths factorize the 2n x 2n Hamiltonian and the sign-function path the 2n x 2n projector
// (I - sign(H)) / 2, so a bound counted at n x n would under-count every path by about a factor of eight.
template <std::size_t NX>
inline constexpr std::int32_t care_extraction_basis_dimension = 2 * static_cast<std::int32_t>(NX);

// 2m^3 for the m x m Householder QR and 3m^3 for applying the orthogonal factor and solving the triangular right-hand
// sides, at m = care_extraction_basis_dimension.
template <std::size_t NX>
inline constexpr std::int32_t care_extraction_rounding_ops = 2 * care_extraction_basis_dimension<NX> * care_extraction_basis_dimension<NX> * care_extraction_basis_dimension<NX> +
        3 * care_extraction_basis_dimension<NX> * care_extraction_basis_dimension<NX> * care_extraction_basis_dimension<NX>;

// On top of the extraction count: 3n^2 for symmetrization, four n x n products at n^2 (2n - 1) each, and three n x n
// sums.
template <std::size_t NX>
constexpr std::int32_t care_residual_rounding_ops()
{
    constexpr std::int32_t n = static_cast<std::int32_t>(NX);
    constexpr std::int32_t extraction_rounding_ops = care_extraction_rounding_ops<NX>;
    constexpr std::int32_t symmetrize_rounding_ops = 3 * n * n;
    constexpr std::int32_t matrix_product_entry_rounding_ops = 2 * n - 1;
    constexpr std::int32_t residual_matrix_products = 4;
    constexpr std::int32_t residual_product_rounding_ops = residual_matrix_products * n * n * matrix_product_entry_rounding_ops;
    constexpr std::int32_t residual_sums = 3;
    constexpr std::int32_t residual_sum_rounding_ops = residual_sums * n * n;
    return extraction_rounding_ops + symmetrize_rounding_ops + residual_product_rounding_ops + residual_sum_rounding_ops;
}

// The scale is the largest term entering the final sum, not the residual itself.
template <typename Scalar, std::size_t NX>
bool care_residual_within_rounding(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NX)>& S, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& P)
{
    using MatN = Eigen::Matrix<Scalar, int(NX), int(NX)>;
    const MatN At_P = (A.transpose() * P).eval();
    const MatN P_A = (P * A).eval();
    const MatN P_S_P = (P * S * P).eval();
    const MatN residual = (At_P + P_A - P_S_P + Q).eval();

    constexpr std::int32_t residual_rounding_ops = care_residual_rounding_ops<NX>();
    const resolved_magnitude<Scalar> residual_scale = largest_magnitude<Scalar>({resolve_magnitude(At_P), resolve_magnitude(P_A), resolve_magnitude(P_S_P), resolve_magnitude(Q)});
    const resolved_magnitude<Scalar> residual_floor = counted_rounding_bound(residual_scale, residual_rounding_ops);
    return magnitude_within(resolve_magnitude(residual), residual_floor);
}

// Only the quasi-triangular factor is read, so the orthogonal factor is not accumulated (about 10n^3 instead of the
// 25n^3 Eigen's default costs).
template <typename Scalar, int N>
bool closed_loop_strictly_left_half_plane(const Eigen::Matrix<Scalar, N, N>& closed_loop)
{
    Eigen::RealSchur<Eigen::Matrix<Scalar, N, N>> closed_loop_schur(closed_loop, false);
    if(closed_loop_schur.info() != Eigen::Success)
        return false;

    const Eigen::Matrix<Scalar, N, N> closed_loop_T = closed_loop_schur.matrixT();
    if(!closed_loop_T.allFinite())
        return false;

    // A backward-stable real Schur factor has an n-operation eigenvalue uncertainty at its largest-entry scale. Its
    // diagonal is not the spectrum: reading it refuses oscillatory loops, and a 2 x 2 block's half-trace can hide a real
    // root in the right half-plane, so each block's root with the larger real part is tested.
    constexpr std::int32_t closed_loop_eigenvalue_rounding_ops = N;
    const resolved_magnitude<Scalar> closed_loop_margin = counted_rounding_bound(resolve_largest_entry(closed_loop_T), closed_loop_eigenvalue_rounding_ops);
    if(!closed_loop_margin.resolved)
        return false;

    return quasi_triangular_spectrum_strictly_left_of<Scalar, N>(closed_loop_T, closed_loop_margin.value);
}

template <typename Scalar, std::size_t NX>
bool care_solution_satisfies_postconditions(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H, const Eigen::Matrix<Scalar, int(NX), int(NX)>& P)
{
    constexpr int n = int(NX);
    using MatN = Eigen::Matrix<Scalar, n, n>;

    const MatN A = H.template block<n, n>(0, 0);
    const MatN S = -H.template block<n, n>(0, n);
    const MatN Q = -H.template block<n, n>(n, 0);
    if(!care_residual_within_rounding<Scalar, NX>(A, S, Q, P))
        return false;

    const MatN closed_loop = (A - S * P).eval();
    return closed_loop_strictly_left_half_plane<Scalar, n>(closed_loop);
}

inline care_error care_error_from_extract(riccati_extract_error error)
{
    switch(error)
    {
        case riccati_extract_error::singular_u11:
            return care_error::singular_u11;
        case riccati_extract_error::non_finite:
            return care_error::non_finite_input;
        case riccati_extract_error::non_psd:
            return care_error::non_psd_solution;
    }
    return care_error::non_finite_input;
}

}

#endif
