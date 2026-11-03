#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_CONDITIONING_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_CONDITIONING_H

#include <Eigen/Core>

#include <concepts>

namespace ctrlpp::detail
{

// Rejects an ill-conditioned swap by |R(last, last)| / |R(0, 0)| of the rank-revealing QR of the Sylvester
// operator, which the swap computes anyway.
struct pivot_ratio_conditioning
{
};

// Hager-Higham 1-norm estimate of the inverse Sylvester operator; not implemented, so selecting it fails to
// compile.
struct hager_higham_conditioning
{
};

// The solvers constrain their trailing policy parameter on this, so overload resolution rejects a foreign
// type in that slot instead of finding an ambiguity with the overloads taking one more matrix.
template <typename T>
concept conditioning_policy = std::same_as<T, pivot_ratio_conditioning> || std::same_as<T, hager_higham_conditioning>;

namespace detail_constants
{

// LAPACK DLANV2 MULTPL: the real-versus-complex-pair decision multiplier.
template <typename Scalar>
constexpr Scalar dlanv2_multpl() noexcept
{
    return Scalar{4};
}

// LAPACK DLAEXC TEN: the post-swap residual threshold multiplier.
template <typename Scalar>
constexpr Scalar rejection_multiplier() noexcept
{
    return Scalar{10};
}

}

// `placed` counts the selected eigenvalues moved to the top; `complete` is false once a swap was rejected;
// `subspace_separation` is the smallest pivot ratio over the swaps, the analogue of LAPACK's SEP.
template <typename Scalar>
struct reorder_result
{
    reorder_result()
        : complete{true}
        , placed{0}
        , subspace_separation{Scalar{1}}
    {
    }

    bool complete;
    Eigen::Index placed;
    Scalar subspace_separation;
};

}

#endif
