#ifndef HPP_GUARD_CTRLPP_DETAIL_HAMILTONIAN_BALANCE_H
#define HPP_GUARD_CTRLPP_DETAIL_HAMILTONIAN_BALANCE_H

// Plain DGEBAL balancing destroys Hamiltonian structure, which is correct here only because the real-Schur
// path treats H as a general 2n x 2n matrix; Benner's symplectic balancing is the structure-preserving
// alternative.
// @cite benner2001    : Benner, "Symplectic Balancing of Hamiltonian Matrices", 2001 (structured alternative, deferred)

#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/diagonal_balance.h"
#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/care_postconditions.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/control/care_types.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
void balance_hamiltonian(Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H, Eigen::Matrix<Scalar, 2 * int(NX), 1>& D_out)
{
    balance_diagonal<Scalar, 2 * int(NX)>(H, D_out);
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

template <typename Scalar, int N, typename Cond>
reorder_result<Scalar> reorder_left_half_plane(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U)
{
    const Scalar lhp_margin = schur_eigenvalue_margin<Scalar, N>(T);
    auto predicate = [lhp_margin](std::complex<Scalar> lam) -> bool { return lam.real() < -lhp_margin; };
    return reorder_real_schur<Scalar, N>(T, U, predicate, Cond{});
}

// Balances a finite copy of H_in into D, then leaves the stable-subspace basis of the balanced Hamiltonian in the
// leading n columns of U.
template <typename Scalar, std::size_t NX, typename Cond>
ctrlpp::expected<reorder_result<Scalar>, care_error> balanced_stable_schur_basis(
        const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in, Eigen::Matrix<Scalar, 2 * int(NX), 1>& D, Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& U)
{
    constexpr int n2 = 2 * int(NX);
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    if(!H_in.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    Mat2N H = H_in;
    balance_hamiltonian<Scalar, NX>(H, D);
    Eigen::RealSchur<Mat2N> schur(H);
    if(schur.info() != Eigen::Success)
        return ctrlpp::unexpected(care_error::schur_failed);

    Mat2N T = schur.matrixT();
    U = schur.matrixU();
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    const reorder_result<Scalar> rr = reorder_left_half_plane<Scalar, n2, Cond>(T, U);
    if(rr.placed < int(NX))
        return ctrlpp::unexpected(care_error::non_lhp_stabilizable);
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);
    return rr;
}

// The left-half-plane invariant subspace of D^-1 H D is D^-1 V, so scaling the leading n columns of U by D
// recovers the caller's basis V.
template <typename Scalar, std::size_t NX, conditioning_policy Cond = pivot_ratio_conditioning>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care_solve_via_balanced_schur(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in)
{
    constexpr int n = int(NX);
    constexpr int n2 = 2 * n;

    Eigen::Matrix<Scalar, n2, 1> D;
    Eigen::Matrix<Scalar, n2, n2> U;
    const ctrlpp::expected<reorder_result<Scalar>, care_error> rr = balanced_stable_schur_basis<Scalar, NX, Cond>(H_in, D, U);
    if(!rr)
        return ctrlpp::unexpected(rr.error());
    U.leftCols(n).array().colwise() *= D.array();

    care_result<Scalar, NX> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if(!P_err)
        return ctrlpp::unexpected(care_error_from_extract(P_err.error()));

    // The back-scaled P solves the caller's equation, not the balanced one, and D equilibrates the norms, so
    // a residual against the balanced H would answer a different question on a different scale.
    if(!care_solution_satisfies_postconditions<Scalar, NX>(H_in, out.P))
        return ctrlpp::unexpected(care_error::unverified_solution);

    out.subspace_separation = rr->subspace_separation;
    out.reorder_complete = rr->complete;
    return out;
}
}

#endif
