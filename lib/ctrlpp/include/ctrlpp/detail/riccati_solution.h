#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H

// Recovers the stabilizing Riccati solution P = U21 U11^-1 in real arithmetic from an orthogonal basis U whose leading
// n columns span the stable invariant subspace.
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/detail/covariance_ops.h"
#include "ctrlpp/detail/psd_acceptance.h"

#include "ctrlpp/expected.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <initializer_list>

namespace ctrlpp::detail
{

// A magnitude whose `value` means something only when `resolved`. Unresolved is contagious and no comparison with an
// unresolved side succeeds, because a plain Frobenius norm overflows to infinity (where inf <= inf certifies anything)
// or underflows to zero (where 0 <= 0 passes): a magnitude the arithmetic could not form declines instead.
template <typename Scalar>
struct resolved_magnitude
{
    Scalar value;
    bool resolved;
};

// The largest-entry rescale max|M| * ||M / max|M|||, which cannot leave the range before its final product. Eigen's
// blueNorm() returns zero for the smallest subnormal, and stableNorm() costs more at the sizes both Riccati paths carry.
template <typename Derived>
resolved_magnitude<typename Derived::Scalar> resolve_magnitude(const Eigen::MatrixBase<Derived>& operand)
{
    using Scalar = typename Derived::Scalar;
    const auto& evaluated = operand.eval();

    if(!evaluated.allFinite())
        return {Scalar{0}, false};

    const Scalar largest_entry = evaluated.cwiseAbs().maxCoeff();
    if(!(largest_entry > Scalar{0}))
        return {Scalar{0}, true};

    const Scalar magnitude = largest_entry * (evaluated / largest_entry).norm();
    if(!std::isfinite(magnitude))
        return {Scalar{0}, false};
    return {magnitude, true};
}

// Squares nothing, so only the finiteness guard is needed. A bound derived from the largest entry differs from the
// Frobenius one, and substituting one for the other would move a threshold.
template <typename Derived>
resolved_magnitude<typename Derived::Scalar> resolve_largest_entry(const Eigen::MatrixBase<Derived>& operand)
{
    using Scalar = typename Derived::Scalar;
    const auto& evaluated = operand.eval();

    if(!evaluated.allFinite())
        return {Scalar{0}, false};
    return {evaluated.cwiseAbs().maxCoeff(), true};
}

template <typename Scalar>
resolved_magnitude<Scalar> largest_magnitude(std::initializer_list<resolved_magnitude<Scalar>> operands)
{
    resolved_magnitude<Scalar> largest{Scalar{0}, true};
    for(const resolved_magnitude<Scalar>& operand : operands)
    {
        if(!operand.resolved)
            return {Scalar{0}, false};
        if(operand.value > largest.value)
            largest.value = operand.value;
    }
    return largest;
}

template <typename Scalar>
resolved_magnitude<Scalar> scaled_magnitude(const resolved_magnitude<Scalar>& operand, Scalar factor)
{
    if(!operand.resolved || !std::isfinite(factor))
        return {Scalar{0}, false};

    const Scalar scaled = operand.value * factor;
    if(!std::isfinite(scaled))
        return {Scalar{0}, false};
    return {scaled, true};
}

template <typename Scalar>
resolved_magnitude<Scalar> counted_rounding_bound(const resolved_magnitude<Scalar>& scale, std::int32_t rounding_ops)
{
    return scaled_magnitude(scale, static_cast<Scalar>(rounding_ops) * std::numeric_limits<Scalar>::epsilon());
}

template <typename Scalar>
bool magnitude_within(const resolved_magnitude<Scalar>& value, const resolved_magnitude<Scalar>& bound)
{
    return value.resolved && bound.resolved && value.value <= bound.value;
}

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

// Relative forward error of a discrete Riccati solution, estimated by inverting the residual map's Frechet derivative.
// At the solution DR_P[E] = A_cl' E A_cl - E with A_cl = A - B K, so with Omega(X) = X - A_cl' X A_cl,
// R(P_hat) = -Omega(E) + O(||E||^2), and the estimate is ||Omega^-1(R(P_hat))||_F / ||P||_F. It is formed on
// P / max|P_ij| and residual / max|P_ij|, under which it is invariant and which cannot leave the range. The residual
// alone cannot decide acceptance: it is a four-term cancellation set by the terms' conditioning, and no threshold on it
// separates accurate answers from inaccurate ones. A first-order estimate under-predicts where the residual cancels
// favorably. A Bartels-Stewart Stein solve is cheaper above N = 4 but let more wrong answers through on the acceptance
// population, so Omega is assembled densely on the symmetric subspace and solved by a rank-revealing QR.
// @cite laub1979  -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 19
template <typename Scalar, int N>
resolved_magnitude<Scalar> estimate_riccati_forward_error(
        const Eigen::Matrix<Scalar, N, N>& closed_loop, const Eigen::Matrix<Scalar, N, N>& residual, const Eigen::Matrix<Scalar, N, N>& P)
{
    static_assert(N > 0, "state dimension must be positive");
    constexpr int M = N * (N + 1) / 2;
    using MatNxN = Eigen::Matrix<Scalar, N, N>;

    if(!closed_loop.allFinite() || !residual.allFinite() || !P.allFinite())
        return {Scalar{0}, false};

    const Scalar peak = P.cwiseAbs().maxCoeff();
    if(!(peak > Scalar{0}))
    {
        // Every residual term containing P is exactly zero, so an exact solve leaves exactly zero.
        const Scalar remainder = residual.cwiseAbs().maxCoeff();
        return {remainder > Scalar{0} ? std::numeric_limits<Scalar>::infinity() : Scalar{0}, true};
    }

    const MatNxN residual_unit = (residual / peak).eval();
    const MatNxN P_unit = (P / peak).eval();
    if(!residual_unit.allFinite())
        return {Scalar{0}, false};

    // Symmetric basis E_ii = e_i e_i', E_ij = e_i e_j' + e_j e_i' for i < j, so Omega(E_ij) = E_ij - c_i c_j' - c_j c_i'
    // with c_i the i-th row of A_cl as a column. The coordinates are the upper triangle in the same order, exact for a
    // symmetric image because the off-diagonal basis element carries the 1 at both places.
    Eigen::Matrix<Scalar, M, M> stein_operator;
    Eigen::Matrix<Scalar, M, 1> stein_rhs;
    Eigen::Index column = 0;
    for(Eigen::Index i = 0; i < N; ++i)
    {
        for(Eigen::Index j = i; j < N; ++j)
        {
            MatNxN image = MatNxN::Zero();
            image(i, j) += Scalar{1};
            const auto ci = closed_loop.row(i).transpose();
            const auto cj = closed_loop.row(j).transpose();
            image -= ci * cj.transpose();
            if(i != j)
            {
                image(j, i) += Scalar{1};
                image -= cj * ci.transpose();
            }

            Eigen::Index row = 0;
            for(Eigen::Index a = 0; a < N; ++a)
                for(Eigen::Index b = a; b < N; ++b)
                    stein_operator(row++, column) = image(a, b);
            ++column;
        }
    }
    Eigen::Index row = 0;
    for(Eigen::Index a = 0; a < N; ++a)
        for(Eigen::Index b = a; b < N; ++b)
        {
            // Only the symmetric part of the residual lies in the operator's range.
            stein_rhs(row++) = (residual_unit(a, b) + residual_unit(b, a)) / Scalar{2};
        }

    // Declared over a Ref to the operator's own storage, the decomposition factorizes in place, so no second M x M array
    // is live; the Householder sequence and the solve are those of the copying form.
    Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::Matrix<Scalar, M, M>>> stein_qr(stein_operator);
    stein_qr.setThreshold(Scalar{M} * std::numeric_limits<Scalar>::epsilon());
    if(!stein_qr.isInvertible())
    {
        // Omega is exactly singular only for a closed-loop pair with lambda_i lambda_j = 1, but the relative rank test also
        // fires once 1 - |lambda|^2 < M eps, where the answer may be good. That is an absence of evidence; the stabilizing
        // claim belongs to the discrete solver's spectrum check.
        return {Scalar{0}, false};
    }

    const Eigen::Matrix<Scalar, M, 1> coordinates = stein_qr.solve(stein_rhs);
    if(!coordinates.allFinite())
        return {Scalar{0}, false};

    MatNxN error_estimate = MatNxN::Zero();
    row = 0;
    for(Eigen::Index a = 0; a < N; ++a)
        for(Eigen::Index b = a; b < N; ++b)
        {
            error_estimate(a, b) = coordinates(row);
            error_estimate(b, a) = coordinates(row);
            ++row;
        }

    const resolved_magnitude<Scalar> numerator = resolve_magnitude(error_estimate);
    const resolved_magnitude<Scalar> denominator = resolve_magnitude(P_unit);
    if(!numerator.resolved || !denominator.resolved || !(denominator.value > Scalar{0}))
        return {Scalar{0}, false};

    const Scalar ratio = numerator.value / denominator.value;
    if(!std::isfinite(ratio))
        return {std::numeric_limits<Scalar>::infinity(), true};
    return {ratio, true};
}

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
