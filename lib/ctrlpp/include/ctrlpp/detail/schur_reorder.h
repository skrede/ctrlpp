#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_REORDER_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_REORDER_H

// Real-Schur block reordering by the Bai and Demmel predicated swap of adjacent diagonal blocks (LAPACK
// DTREXC, DLAEXC), with DLANV2 Murnaghan standardization of every 2x2 block a swap touches.
// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/detail/schur_swap.h"
#include "ctrlpp/detail/quasi_triangular.h"
#include "ctrlpp/detail/schur_standardize.h"
#include "ctrlpp/detail/schur_conditioning.h"

#include <Eigen/Dense>

#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

// `first` is the larger real root of a real pair and shares its modulus and real part with the other root of
// a conjugate pair, so the predicate asks about `first` alone.
template <typename Scalar, int N, typename Predicate>
bool schur_block_matches(const Eigen::Matrix<Scalar, N, N>& T, Eigen::Index pos, Eigen::Index n_block, Predicate& predicate)
{
    return predicate(quasi_triangular_block_spectrum<Scalar, N>(T, pos, n_block).first);
}

template <typename Scalar, int N, typename Predicate>
bool find_matching_block(const Eigen::Matrix<Scalar, N, N>& T, Scalar scale_T, Eigen::Index& scan, Eigen::Index& scan_nb, Predicate& predicate)
{
    while(scan < N)
    {
        const Eigen::Index nb = quasi_triangular_block_size<Scalar, N>(T, scan, scale_T);
        if(schur_block_matches<Scalar, N>(T, scan, nb, predicate))
        {
            scan_nb = nb;
            return true;
        }
        scan += nb;
    }
    return false;
}

// A rejected swap leaves the block where it is and marks the result incomplete.
template <typename Scalar, int N, typename Cond>
bool swap_with_left_block(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Scalar scale_T, Eigen::Index& cur, Eigen::Index scan_nb, Cond tag, reorder_result<Scalar>& r)
{
    // A 2x2 block ending at cur - 1 starts at cur - 2.
    const Eigen::Index left_nb = (cur >= 2) ? quasi_triangular_block_size<Scalar, N>(T, cur - 2, scale_T) : 1;
    const Eigen::Index left_pos = cur - left_nb;

    Scalar pivot_ratio = Scalar{1};
    const bool ok = swap_real_schur_blocks<Scalar, N, Cond>(T, U, left_pos, left_nb, scan_nb, tag, pivot_ratio);
    if(!ok)
    {
        r.complete = false;
        r.subspace_separation = std::min(r.subspace_separation, pivot_ratio);
        return false;
    }
    r.subspace_separation = std::min(r.subspace_separation, pivot_ratio);

    if(scan_nb == 2)
        standardize_2x2_block<Scalar, N>(T, U, left_pos);
    if(left_nb == 2)
        standardize_2x2_block<Scalar, N>(T, U, left_pos + scan_nb);
    cur = left_pos;
    return true;
}

// Standardization can split a 2x2 block whose conjugate pair degenerated into two real eigenvalues, and the
// split block may no longer match, which stops the bubble.
template <typename Scalar, int N, typename Predicate, typename Cond>
Eigen::Index bubble_block_up(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Scalar scale_T, Eigen::Index placed_pos, Eigen::Index cur, Eigen::Index& scan_nb,
        Predicate& predicate, Cond tag, reorder_result<Scalar>& r)
{
    while(cur > placed_pos)
    {
        if(!swap_with_left_block<Scalar, N, Cond>(T, U, scale_T, cur, scan_nb, tag, r))
            break;
        scan_nb = quasi_triangular_block_size<Scalar, N>(T, cur, scale_T);
        if(!schur_block_matches<Scalar, N>(T, cur, scan_nb, predicate))
            break;
    }
    return cur;
}

// The block sizes are judged against the scale of the factor as it arrived. A bubble that stops short resumes
// the scan past the stuck block, which guarantees termination.
template <typename Scalar, int N, typename Predicate, typename Cond = pivot_ratio_conditioning>
reorder_result<Scalar> reorder_real_schur(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Predicate&& predicate, Cond tag = {})
{
    static_assert(std::is_same_v<Cond, pivot_ratio_conditioning>, "Cond must be pivot_ratio_conditioning; hager_higham_conditioning is not implemented");

    reorder_result<Scalar> r;
    const Scalar scale_T = T.cwiseAbs().maxCoeff();
    Eigen::Index placed_pos = 0;
    while(placed_pos < N)
    {
        Eigen::Index scan = placed_pos;
        Eigen::Index scan_nb = 1;
        if(!find_matching_block<Scalar, N>(T, scale_T, scan, scan_nb, predicate))
            break;

        const Eigen::Index cur = bubble_block_up<Scalar, N>(T, U, scale_T, placed_pos, scan, scan_nb, predicate, tag, r);
        if(cur == placed_pos)
        {
            r.placed += scan_nb;
            placed_pos += scan_nb;
        }
        else
            placed_pos = cur + quasi_triangular_block_size<Scalar, N>(T, cur, scale_T);
    }
    return r;
}

}

#endif
