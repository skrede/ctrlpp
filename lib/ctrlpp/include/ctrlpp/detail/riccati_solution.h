#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H

// Recovers the stabilizing Riccati solution P = U21 U11^-1 in real arithmetic from an orthogonal basis U whose leading
// n columns span the stable invariant subspace.
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/detail/covariance_ops.h"
#include "ctrlpp/detail/psd_acceptance.h"
#include "ctrlpp/detail/riccati_magnitude.h"
#include "ctrlpp/detail/riccati_forward_error.h"

#include "ctrlpp/expected.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

// Rounded operations along the longest chain producing one entry of the discrete Riccati residual: six contractions over
// the state dimension (A'P, (A'P)A, B'P, (B'P)A, A'PB, and A'PB against the gain) at 2 NX - 1 each, two over the input
// dimension ((B'P)B and the gain's inner dimension) at 2 NU - 1 each, one addition for R + B'PB, 2 NU for the gain's
// rank-revealing QR solve, and three additions assembling the terms. Every operation is counted whether or not it rounds.
template <std::size_t NX, std::size_t NU>
constexpr std::int32_t dare_residual_ops = 6 * (2 * static_cast<std::int32_t>(NX) - 1) + 2 * (2 * static_cast<std::int32_t>(NU) - 1) + 1 + 2 * static_cast<std::int32_t>(NU) + 3;

// Two posings that differ only by a positive common weight scale share one gain exactly (the gain is homogeneous of
// degree zero in (P, Q, R)), so only rounding separates them. The square root is taken because two independently
// rounded results are compared with each other rather than with the exact value.
template <typename Scalar, std::size_t NX, std::size_t NU>
Scalar dare_gain_agreement_margin()
{
    return std::sqrt(Scalar{dare_residual_ops<NX, NU>} * std::numeric_limits<Scalar>::epsilon());
}

// For a radix-2 type with eps = 2^-p, a relative accuracy of sqrt(eps) = 2^(-p/2) retains exactly p/2 of the p
// significand bits, so the margin follows from the type rather than from a calibration.
template <typename Scalar>
Scalar half_significand_margin()
{
    return std::sqrt(std::numeric_limits<Scalar>::epsilon());
}

// An unresolved estimate has produced no evidence, which is neither a certificate nor a refutation.
enum class riccati_accuracy
{
    within,
    exceeded,
    unresolved,
};

template <typename Scalar, int N>
riccati_accuracy riccati_forward_error_verdict(const Eigen::Matrix<Scalar, N, N>& closed_loop, const Eigen::Matrix<Scalar, N, N>& residual, const Eigen::Matrix<Scalar, N, N>& P)
{
    const resolved_magnitude<Scalar> estimate = estimate_riccati_forward_error<Scalar, N>(closed_loop, residual, P);
    if(!estimate.resolved)
        return riccati_accuracy::unresolved;
    return estimate.value <= half_significand_margin<Scalar>() ? riccati_accuracy::within : riccati_accuracy::exceeded;
}

enum class riccati_extract_error
{
    singular_u11,
    non_finite,
    non_psd,
};

// Writes into P_out so a hot-path caller avoids copying the matrix out of an expected return.
template <typename Scalar, int N2>
ctrlpp::expected<void, riccati_extract_error> extract_riccati_solution_into(Eigen::Matrix<Scalar, N2 / 2, N2 / 2>& P_out, const Eigen::Matrix<Scalar, N2, N2>& U)
{
    static_assert(N2 > 0 && (N2 % 2 == 0), "U must have even size 2n x 2n");

    constexpr int n = N2 / 2;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;

    const MatNxN U11 = U.template block<n, n>(0, 0);
    const MatNxN U21 = U.template block<n, n>(n, 0);

    auto qr_U11T = U11.transpose().colPivHouseholderQr();
    if(!qr_U11T.isInvertible())
        return ctrlpp::unexpected(riccati_extract_error::singular_u11);

    // P = U21 U11^-1 <=> U11^T P^T = U21^T, solved against the QR of U11^T rather than through an inverse.
    const MatNxN P_raw = qr_U11T.solve(U21.transpose()).transpose();
    P_out = ctrlpp::detail::symmetrize(P_raw);

    if(!P_out.allFinite())
        return ctrlpp::unexpected(riccati_extract_error::non_finite);

    if(!psd_within_rounding<Scalar, n>(P_out))
        return ctrlpp::unexpected(riccati_extract_error::non_psd);

    return {};
}

template <typename Scalar, int N2>
ctrlpp::expected<Eigen::Matrix<Scalar, N2 / 2, N2 / 2>, riccati_extract_error> extract_riccati_solution(const Eigen::Matrix<Scalar, N2, N2>& U)
{
    Eigen::Matrix<Scalar, N2 / 2, N2 / 2> P;
    auto err = extract_riccati_solution_into<Scalar, N2>(P, U);
    if(!err)
        return ctrlpp::unexpected(err.error());
    return P;
}

}

#endif
