#ifndef HPP_GUARD_CTRLPP_DETAIL_PSD_ACCEPTANCE_H
#define HPP_GUARD_CTRLPP_DETAIL_PSD_ACCEPTANCE_H

#include "ctrlpp/detail/power_of_two_scale.h"

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include <limits>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

// The diagonal shift tau that lets a Cholesky factorization complete on every symmetric matrix within one rounding per
// entry of a positive semi-definite one, given its Frobenius norm f as `norm()` evaluates it and d = max(0, max_i p_ii).
// tau exceeds [u f / (1 - u) + (u + c_n (1 + u)) d] / [1 - u - c_n (1 + u)]: one rounding per entry moves the smallest
// eigenvalue by at most u f / (1 - u) (Weyl), adding tau rounds each diagonal entry once, and the completion condition
// is c_n = n gamma_{n+1} / (1 - gamma_{n+1}), gamma_k = k u / (1 - k u). The evaluation, the norm's included, carries
// m = ceil(n^2 / 2) + 7 roundings, so dividing by 1 - (m + 1) u leaves tau strictly above that bound.
// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Thm 10.7 (p. 200), Lemma 3.1
// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1 (Weyl)
template <typename Scalar, int N>
Scalar cholesky_completion_shift(Scalar frobenius, Scalar largest_diagonal)
{
    constexpr Scalar u = std::numeric_limits<Scalar>::epsilon() / 2;
    constexpr Scalar n = static_cast<Scalar>(N);
    static_assert(N > 0 && n * (n + 1) * u <= Scalar{1} / 8, "the shift's rounding count assumes a fixed size with n (n + 1) u <= 1/8");
    constexpr std::int32_t roundings = (N * N + 1) / 2 + 7;

    // n gamma_{n+1} / (1 - gamma_{n+1}) = n (n + 1) u / (1 - 2 (n + 1) u), whose operands are exact.
    const Scalar completion = n * (n + 1) * u / (1 - 2 * (n + 1) * u);
    const Scalar weight = u + completion + completion * u;
    const Scalar shift = (frobenius * (u / (1 - u)) + weight * largest_diagonal) / (1 - weight);
    return shift / (1 - static_cast<Scalar>(roundings + 1) * u);
}

// Whether a symmetric matrix is within one rounding per entry of a positive semi-definite matrix. The check runs on P
// divided by the power of two at or above its largest entry, so the verdict is invariant under exact power-of-two
// scaling and nothing it forms can overflow. There a Cholesky factorization of P + tau I completes on every such matrix
// (Demmel's completion condition) and fails on any whose smallest eigenvalue is below -(tau + (c_n (1 + u) + u) (d + tau)).
// Only the lower triangle is read.
// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Thm 10.3 (p. 197), Thm 10.7 (p. 200)
template <typename Scalar, int N>
bool psd_within_rounding(const Eigen::Matrix<Scalar, N, N>& P)
{
    if(!P.allFinite())
        return false;
    const Scalar largest = P.cwiseAbs().maxCoeff();
    if(largest == Scalar{0})
        return true;

    Eigen::Matrix<Scalar, N, N> shifted = P / power_of_two_ceiling(largest);
    const Scalar largest_diagonal = std::max(Scalar{0}, shifted.diagonal().maxCoeff());
    shifted.diagonal().array() += cholesky_completion_shift<Scalar, N>(shifted.norm(), largest_diagonal);
    return Eigen::LLT<Eigen::Matrix<Scalar, N, N>>(shifted).info() == Eigen::Success;
}

}

#endif
