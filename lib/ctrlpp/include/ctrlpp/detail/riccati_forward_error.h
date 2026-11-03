#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_FORWARD_ERROR_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_FORWARD_ERROR_H

// Relative forward error of a discrete Riccati solution, estimated by inverting the residual map's Frechet derivative.
// At the solution DR_P[E] = A_cl' E A_cl - E with A_cl = A - B K, so with Omega(X) = X - A_cl' X A_cl,
// R(P_hat) = -Omega(E) + O(||E||^2), and the estimate is ||Omega^-1(R(P_hat))||_F / ||P||_F. It is formed on
// P / max|P_ij| and residual / max|P_ij|, under which it is invariant and which cannot leave the range. The residual
// alone cannot decide acceptance: it is a four-term cancellation set by the terms' conditioning, and no threshold on it
// separates accurate answers from inaccurate ones. A first-order estimate under-predicts where the residual cancels
// favorably. A Bartels-Stewart Stein solve is cheaper at N = 6 and 8 but let more wrong answers through on the acceptance
// population, so Omega is assembled densely on the symmetric subspace and solved by a rank-revealing QR.
// @cite laub1979  -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
// @cite higham2008 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2008, Ch. 19

#include "ctrlpp/detail/riccati_magnitude.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>

namespace ctrlpp::detail
{

template <int N>
inline constexpr int symmetric_dimension = N * (N + 1) / 2;

// Symmetric basis E_ii = e_i e_i', E_ij = e_i e_j' + e_j e_i' for i < j, so Omega(E_ij) = E_ij - c_i c_j' - c_j c_i' with
// c_i the i-th row of A_cl as a column.
template <typename Scalar, int N>
Eigen::Matrix<Scalar, N, N> stein_basis_image(const Eigen::Matrix<Scalar, N, N>& closed_loop, Eigen::Index i, Eigen::Index j)
{
    using MatNxN = Eigen::Matrix<Scalar, N, N>;
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
    return image;
}

// The coordinates are the upper triangle in basis order, exact for a symmetric image because the off-diagonal basis
// element carries the 1 at both places.
template <typename Scalar, int N>
void assemble_stein_operator(const Eigen::Matrix<Scalar, N, N>& closed_loop, Eigen::Matrix<Scalar, symmetric_dimension<N>, symmetric_dimension<N>>& stein_operator)
{
    Eigen::Index column = 0;
    for(Eigen::Index i = 0; i < N; ++i)
    {
        for(Eigen::Index j = i; j < N; ++j)
        {
            const Eigen::Matrix<Scalar, N, N> image = stein_basis_image<Scalar, N>(closed_loop, i, j);
            Eigen::Index row = 0;
            for(Eigen::Index a = 0; a < N; ++a)
                for(Eigen::Index b = a; b < N; ++b)
                    stein_operator(row++, column) = image(a, b);
            ++column;
        }
    }
}

template <typename Scalar, int N>
Eigen::Matrix<Scalar, symmetric_dimension<N>, 1> stein_right_hand_side(const Eigen::Matrix<Scalar, N, N>& residual_unit)
{
    Eigen::Matrix<Scalar, symmetric_dimension<N>, 1> stein_rhs;
    Eigen::Index row = 0;
    for(Eigen::Index a = 0; a < N; ++a)
        for(Eigen::Index b = a; b < N; ++b)
        {
            // Only the symmetric part of the residual lies in the operator's range.
            stein_rhs(row++) = (residual_unit(a, b) + residual_unit(b, a)) / Scalar{2};
        }
    return stein_rhs;
}

// Declared over a Ref to the operator's own storage, the decomposition factorizes in place, so no second M x M array is
// live; the Householder sequence and the solve are those of the copying form.
template <typename Scalar, int M>
bool solve_stein_system(Eigen::Matrix<Scalar, M, M>& stein_operator, const Eigen::Matrix<Scalar, M, 1>& stein_rhs, Eigen::Matrix<Scalar, M, 1>& coordinates)
{
    Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::Matrix<Scalar, M, M>>> stein_qr(stein_operator);
    stein_qr.setThreshold(Scalar{M} * std::numeric_limits<Scalar>::epsilon());
    // Omega is exactly singular only for a closed-loop pair with lambda_i lambda_j = 1, but the relative rank test also
    // fires once 1 - |lambda|^2 < M eps, where the answer may be good. That is an absence of evidence; the stabilizing
    // claim belongs to the discrete solver's spectrum check.
    if(!stein_qr.isInvertible())
        return false;

    coordinates = stein_qr.solve(stein_rhs);
    return coordinates.allFinite();
}

template <typename Scalar, int N>
Eigen::Matrix<Scalar, N, N> unpack_symmetric(const Eigen::Matrix<Scalar, symmetric_dimension<N>, 1>& coordinates)
{
    Eigen::Matrix<Scalar, N, N> unpacked = Eigen::Matrix<Scalar, N, N>::Zero();
    Eigen::Index row = 0;
    for(Eigen::Index a = 0; a < N; ++a)
        for(Eigen::Index b = a; b < N; ++b)
        {
            unpacked(a, b) = coordinates(row);
            unpacked(b, a) = coordinates(row);
            ++row;
        }
    return unpacked;
}

template <typename Scalar, int N>
resolved_magnitude<Scalar> forward_error_ratio(const Eigen::Matrix<Scalar, N, N>& error_estimate, const Eigen::Matrix<Scalar, N, N>& P_unit)
{
    const resolved_magnitude<Scalar> numerator = resolve_magnitude(error_estimate);
    const resolved_magnitude<Scalar> denominator = resolve_magnitude(P_unit);
    if(!numerator.resolved || !denominator.resolved || !(denominator.value > Scalar{0}))
        return {Scalar{0}, false};

    const Scalar ratio = numerator.value / denominator.value;
    if(!std::isfinite(ratio))
        return {std::numeric_limits<Scalar>::infinity(), true};
    return {ratio, true};
}

// With P exactly zero every residual term containing P is too, so an exact solve leaves exactly zero.
template <typename Scalar, int N>
resolved_magnitude<Scalar> zero_solution_forward_error(const Eigen::Matrix<Scalar, N, N>& residual)
{
    const Scalar remainder = residual.cwiseAbs().maxCoeff();
    return {remainder > Scalar{0} ? std::numeric_limits<Scalar>::infinity() : Scalar{0}, true};
}

template <typename Scalar, int N>
resolved_magnitude<Scalar> estimate_riccati_forward_error(
        const Eigen::Matrix<Scalar, N, N>& closed_loop, const Eigen::Matrix<Scalar, N, N>& residual, const Eigen::Matrix<Scalar, N, N>& P)
{
    static_assert(N > 0, "state dimension must be positive");
    constexpr int M = symmetric_dimension<N>;
    using MatNxN = Eigen::Matrix<Scalar, N, N>;

    if(!closed_loop.allFinite() || !residual.allFinite() || !P.allFinite())
        return {Scalar{0}, false};

    const Scalar peak = P.cwiseAbs().maxCoeff();
    if(!(peak > Scalar{0}))
        return zero_solution_forward_error<Scalar, N>(residual);

    const MatNxN residual_unit = (residual / peak).eval();
    const MatNxN P_unit = (P / peak).eval();
    if(!residual_unit.allFinite())
        return {Scalar{0}, false};

    Eigen::Matrix<Scalar, M, M> stein_operator;
    assemble_stein_operator<Scalar, N>(closed_loop, stein_operator);
    Eigen::Matrix<Scalar, M, 1> coordinates;
    if(!solve_stein_system<Scalar, M>(stein_operator, stein_right_hand_side<Scalar, N>(residual_unit), coordinates))
        return {Scalar{0}, false};
    return forward_error_ratio<Scalar, N>(unpack_symmetric<Scalar, N>(coordinates), P_unit);
}

}

#endif
