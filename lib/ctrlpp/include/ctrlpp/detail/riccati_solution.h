#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H

/// @brief Riccati invariant-subspace P extraction: P = U21 * U11^-1 on real U.
///
/// Given a reordered orthogonal real basis U of size 2n x 2n (the output of
/// ctrlpp::detail::reorder_real_schur), this header provides the single
/// routine required to recover the stabilising Riccati solution P.
/// P is computed as U21 * U11^-1 directly on real arithmetic, symmetrised
/// via ctrlpp::detail::symmetrize, and checked for finiteness and positive
/// semi-definiteness. The positive semi-definiteness floor is derived from
/// std::numeric_limits<Scalar>::epsilon() scaled by the infinity norm of P;
/// no hardcoded numerical literal appears anywhere in the primitive.
///
/// Errors surface through ctrlpp::expected with a three-variant enum shared
/// between DARE and CARE; the public dare_error / care_error enums map
/// onto these variants one-to-one in the control/ wrappers.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

namespace ctrlpp::detail
{

/// @brief Error variants produced by extract_riccati_solution.
///
/// Shared between DARE and CARE. The public dare_error / care_error enums
/// map their singular_u11 / non_finite_input / non_psd_solution variants
/// directly onto these values.
enum class riccati_extract_error
{
    singular_u11,
    non_finite,
    non_psd,
};

/// @brief Extract the stabilising Riccati solution into a caller-supplied matrix.
///
/// Given an orthogonal real U of size 2n x 2n with the stable invariant
/// subspace in its leading n columns (output of reorder_real_schur), compute
/// P = U21 * U11^-1, symmetrise in place, and validate positive semi-definiteness.
/// Writes directly into P_out so that callers can avoid the ctrlpp::expected<Matrix>
/// return-by-value copy on the hot path.
///
/// The implementation prefers a back-substitution against U11 (solving
/// U11^T * P^T = U21^T) over an explicit inverse for numerical accuracy.
///
/// The positive semi-definiteness floor is -eps * ||P||_inf, matching the
/// LAPACK convention of scaling relative thresholds by the operand norm.
///
/// @returns ctrlpp::expected<void, riccati_extract_error>.
template <typename Scalar, int N2>
[[nodiscard]] auto extract_riccati_solution_into(
    Eigen::Matrix<Scalar, N2 / 2, N2 / 2>&    P_out,
    const Eigen::Matrix<Scalar, N2, N2>&      U)
    -> ctrlpp::expected<void, riccati_extract_error>
{
    static_assert(N2 > 0 && (N2 % 2 == 0), "U must have even size 2n x 2n");

    constexpr int n = N2 / 2;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;

    const MatNxN U11 = U.template block<n, n>(0, 0);
    const MatNxN U21 = U.template block<n, n>(n, 0);

    auto qr_U11T = U11.transpose().colPivHouseholderQr();
    if (!qr_U11T.isInvertible())
        return ctrlpp::unexpected(riccati_extract_error::singular_u11);

    // P = U21 * U11^-1  <=>  U11^T * P^T = U21^T  (solve against QR of U11^T).
    const MatNxN P_raw = qr_U11T.solve(U21.transpose()).transpose();
    P_out = ctrlpp::detail::symmetrize(P_raw);

    if (!P_out.allFinite())
        return ctrlpp::unexpected(riccati_extract_error::non_finite);

    // PSD check via LDLT: by Sylvester's law of inertia, the signs of the pivots
    // in the D diagonal match the signs of the eigenvalues of a symmetric P.
    // Detects non-PSD at ~1/10 the instruction cost of a full eigendecomposition
    // while preserving the same eps-scaled rejection threshold.
    Eigen::LDLT<MatNxN> ldlt(P_out);
    const Scalar psd_floor =
        -std::numeric_limits<Scalar>::epsilon() * P_out.cwiseAbs().maxCoeff();
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() < psd_floor)
        return ctrlpp::unexpected(riccati_extract_error::non_psd);

    return {};
}

/// @brief Value-returning wrapper around `extract_riccati_solution_into`.
template <typename Scalar, int N2>
[[nodiscard]] auto extract_riccati_solution(const Eigen::Matrix<Scalar, N2, N2>& U)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, N2 / 2, N2 / 2>, riccati_extract_error>
{
    Eigen::Matrix<Scalar, N2 / 2, N2 / 2> P;
    auto err = extract_riccati_solution_into<Scalar, N2>(P, U);
    if (!err)
        return ctrlpp::unexpected(err.error());
    return P;
}

}

#endif
