#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_SWAP_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_SWAP_H

// Bai and Demmel's swap of adjacent diagonal blocks of a real Schur factor (LAPACK DLAEXC).
// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include "ctrlpp/detail/schur_standardize.h"
#include "ctrlpp/detail/schur_swap_stages.h"
#include "ctrlpp/detail/schur_conditioning.h"

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <algorithm>

namespace ctrlpp::detail
{

// The 1x1/1x1 swap always succeeds (Bai and Demmel, p. 79). Coinciding eigenvalues make it a no-op, as in
// LAPACK; otherwise (cs, sn) = (t12, -diff) / hypot(t12, diff) zeroes the new subdiagonal without overflow
// (Golub and Van Loan, Alg. 5.1.3).
template <typename Scalar, int N>
bool swap_real_schur_1x1(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Scalar& pivot_ratio_out)
{
    using std::abs;
    using std::hypot;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();

    const Scalar t11 = T(p, p);
    const Scalar t22 = T(p + 1, p + 1);
    const Scalar t12 = T(p, p + 1);
    const Scalar diff = t11 - t22;

    if(abs(diff) < eps * std::max(abs(t11), abs(t22)))
    {
        pivot_ratio_out = Scalar{1};
        return true;
    }

    const Scalar r = hypot(t12, diff);
    const Scalar cs = t12 / r;
    const Scalar sn = -diff / r;
    apply_givens_similarity<Scalar, N>(T, U, p, cs, sn);
    T(p + 1, p) = Scalar{0};
    pivot_ratio_out = Scalar{1};
    return true;
}

// The 1x2, 2x1 and 2x2 swaps share one path through a fixed 4x4 Sylvester operator. The swap is rejected on
// the operator's conditioning or on the post-swap residual (Bai and Demmel, p. 79).
template <typename Scalar, int N>
bool swap_real_schur_2x2_general(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Eigen::Index n1, Eigen::Index n2, Scalar& pivot_ratio_out)
{
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Eigen::Index w = n1 + n2;

    Eigen::Matrix<Scalar, 4, 4> K = Eigen::Matrix<Scalar, 4, 4>::Zero();
    Eigen::Matrix<Scalar, 4, 1> rhs = Eigen::Matrix<Scalar, 4, 1>::Zero();
    assemble_sylvester_operator<Scalar, N>(T, p, n1, n2, K, rhs);

    Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, 4, 4>> qr(K);
    pivot_ratio_out = sylvester_pivot_ratio<Scalar>(qr, n1 * n2);
    if(pivot_ratio_out < eps)
        return false;

    Eigen::Matrix<Scalar, 4, 4> Q = Eigen::Matrix<Scalar, 4, 4>::Identity();
    build_swap_basis<Scalar>(qr.solve(rhs), n1, n2, Q);
    Eigen::Matrix<Scalar, 4, 4> A_new = Eigen::Matrix<Scalar, 4, 4>::Zero();
    tentative_swap_block<Scalar, N>(T, p, w, Q, A_new);
    const Scalar reject_threshold = detail_constants::rejection_multiplier<Scalar>() * eps * T.cwiseAbs().maxCoeff();
    if(swap_residual_norm<Scalar>(A_new, n2, w) > reject_threshold)
        return false;

    write_swapped_block<Scalar, N>(T, p, n2, w, A_new);
    apply_swap_outside_block<Scalar, N>(T, p, w, Q);
    accumulate_swap_basis<Scalar, N>(U, p, w, Q);
    return true;
}

template <typename Scalar, int N, typename Cond>
bool swap_real_schur_blocks(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Eigen::Index n1, Eigen::Index n2, Cond, Scalar& pivot_ratio_out)
{
    if(n1 == 1 && n2 == 1)
        return swap_real_schur_1x1<Scalar, N>(T, U, p, pivot_ratio_out);
    return swap_real_schur_2x2_general<Scalar, N>(T, U, p, n1, n2, pivot_ratio_out);
}

}

#endif
