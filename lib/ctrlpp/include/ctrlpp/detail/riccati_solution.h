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
/// Errors surface through std::expected with a three-variant enum shared
/// between DARE and CARE; the public dare_error / care_error enums map
/// onto these variants one-to-one in the control/ wrappers.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <expected>

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

/// @brief Extract the stabilising Riccati solution from a reordered invariant subspace.
///
/// Given an orthogonal real U of size 2n x 2n with the stable invariant
/// subspace in its leading n columns (output of reorder_real_schur), compute
/// P = U21 * U11^-1, symmetrise, and validate positive semi-definiteness.
///
/// The implementation prefers a back-substitution against U11 (solving
/// U11^T * P^T = U21^T) over an explicit inverse for numerical accuracy:
/// this is a single pass against the column-pivoted QR factor instead of
/// forming and multiplying U11^-1.
///
/// The positive semi-definiteness floor is -eps * ||P||_inf, matching the
/// LAPACK convention of scaling relative thresholds by the operand norm.
///
/// @returns P on success; std::unexpected(riccati_extract_error) otherwise.
template <typename Scalar, int N2>
[[nodiscard]] auto extract_riccati_solution(const Eigen::Matrix<Scalar, N2, N2>& U)
    -> std::expected<Eigen::Matrix<Scalar, N2 / 2, N2 / 2>, riccati_extract_error>
{
    static_assert(N2 > 0 && (N2 % 2 == 0), "U must have even size 2n x 2n");

    constexpr int n = N2 / 2;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;

    const MatNxN U11 = U.template block<n, n>(0, 0);
    const MatNxN U21 = U.template block<n, n>(n, 0);

    auto qr_U11T = U11.transpose().colPivHouseholderQr();
    if (!qr_U11T.isInvertible())
        return std::unexpected(riccati_extract_error::singular_u11);

    // Prefer solve over explicit inverse for numerical accuracy.
    // P = U21 * U11^-1  <=>  U11^T * P^T = U21^T  (solve against QR of U11^T).
    const MatNxN P_raw = qr_U11T.solve(U21.transpose()).transpose();
    const MatNxN P = ctrlpp::detail::symmetrize(P_raw);

    if (!P.allFinite())
        return std::unexpected(riccati_extract_error::non_finite);

    Eigen::SelfAdjointEigenSolver<MatNxN> eigs(P, Eigen::EigenvaluesOnly);
    const Scalar psd_floor =
        -std::numeric_limits<Scalar>::epsilon() * P.cwiseAbs().maxCoeff();
    for (int i = 0; i < n; ++i)
    {
        if (eigs.eigenvalues()(i) < psd_floor)
            return std::unexpected(riccati_extract_error::non_psd);
    }

    return P;
}

}

#endif
