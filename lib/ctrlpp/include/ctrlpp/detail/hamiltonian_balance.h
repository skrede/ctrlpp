#ifndef HPP_GUARD_CTRLPP_DETAIL_HAMILTONIAN_BALANCE_H
#define HPP_GUARD_CTRLPP_DETAIL_HAMILTONIAN_BALANCE_H

/// @brief DGEBAL-style diagonal balancing and balanced-Schur CARE solve.
///
/// `balance_hamiltonian` iteratively equilibrates row and column infinity
/// norms of the 2n x 2n Hamiltonian by applying a diagonal similarity
/// H' = D^{-1} H D, mirroring LAPACK DGEBAL phase 2 on a fixed-size
/// Eigen::Matrix<Scalar, 2*NX, 2*NX> without heap allocation. The LAPACK
/// permutation/isolation phase is omitted because on the Hamiltonians
/// produced by `build_care_hamiltonian` it rarely fires and adds branches
/// without measurable gain. Base-2 scaling constants SCLFAC = 2 and
/// FACTOR = 19/20 (= 0.95) are expressed as rational literals so no bare
/// non-structural numeric constant appears in the hot path.
///
/// `care_solve_via_balanced_schur` composes the balance step, the shared
/// `Eigen::RealSchur` + `ctrlpp::detail::reorder_real_schur` pipeline,
/// the D back-scale on the leading n columns of U, and the shared
/// `ctrlpp::detail::extract_riccati_solution_into` primitive. The back-
/// scale is a left-multiplication by diag(D) on the stable-subspace basis
/// columns, equivalent to `U.leftCols(n).array().colwise() *= D.array()`.
/// This is required because the LHP invariant subspace of H' = D^{-1} H D
/// is D^{-1} V (where V is the LHP invariant subspace of H), so recovering
/// V costs one element-wise product.
///
/// Plain DGEBAL destroys Hamiltonian structure; the structure-preserving
/// Benner 2001 symplectic balancing is reserved for a structured-Schur
/// path (conditional, not wired here). For the plain-Schur pipeline the
/// downstream treats H as a general 2n x 2n matrix, so plain DGEBAL is
/// correct.
///
/// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 algorithm
/// @cite benner2001    : Benner, "Symplectic Balancing of Hamiltonian Matrices", 2001 (structured alternative, deferred)

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <expected>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

/// @brief In-place DGEBAL-style diagonal balance of a 2n x 2n matrix.
///
/// Overwrites H with the balanced matrix H' = D^{-1} H D and writes the
/// diagonal scaling vector into D_out (with D_out(i) >= 0 for all i).
template <typename Scalar, std::size_t NX>
auto balance_hamiltonian(
    Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H,
    Eigen::Matrix<Scalar, 2 * int(NX), 1>&           D_out)
    -> void
{
    constexpr int n2 = 2 * int(NX);
    using Vec2N = Eigen::Matrix<Scalar, n2, 1>;

    const Scalar sclfac     = Scalar{2};
    const Scalar factor_num = Scalar{19};
    const Scalar factor_den = Scalar{20};
    const Scalar factor     = factor_num / factor_den;

    D_out = Vec2N::Ones();
    bool noconv = true;

    while (noconv)
    {
        noconv = false;
        for (int i = 0; i < n2; ++i)
        {
            Scalar r = Scalar{0};
            Scalar c = Scalar{0};
            for (int j = 0; j < n2; ++j)
            {
                if (j == i) continue;
                r += std::abs(H(i, j));
                c += std::abs(H(j, i));
            }
            if (r == Scalar{0} || c == Scalar{0}) continue;

            Scalar g = r / sclfac;
            Scalar f = Scalar{1};
            Scalar s = c + r;

            while (c < g)
            {
                f *= sclfac;
                c *= sclfac * sclfac;
            }
            g = r * sclfac;
            while (c >= g)
            {
                f /= sclfac;
                c /= sclfac * sclfac;
            }

            if ((c + r) < factor * s)
            {
                noconv = true;
                D_out(i) *= f;
                H.row(i) /= f;
                H.col(i) *= f;
            }
        }
    }
}

/// @brief CARE solve via DGEBAL pre-balance then real Schur + Bai-Demmel reorder.
template <typename Scalar, std::size_t NX,
          conditioning_policy Cond = pivot_ratio_conditioning>
auto care_solve_via_balanced_schur(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in)
    -> std::expected<care_result<Scalar, NX>, care_error>
{
    constexpr int n  = int(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;
    using Vec2N = Eigen::Matrix<Scalar, n2, 1>;

    if (!H_in.allFinite())
        return std::unexpected(care_error::non_finite_input);

    Mat2N H = H_in;
    Vec2N D;
    balance_hamiltonian<Scalar, NX>(H, D);

    Eigen::RealSchur<Mat2N> schur(H);
    if (schur.info() != Eigen::Success)
        return std::unexpected(care_error::schur_failed);

    Mat2N T = schur.matrixT();
    Mat2N U = schur.matrixU();
    if (!T.allFinite() || !U.allFinite())
        return std::unexpected(care_error::non_finite_input);

    const Scalar scale      = T.cwiseAbs().maxCoeff();
    const Scalar eps        = std::numeric_limits<Scalar>::epsilon();
    const Scalar lhp_margin = eps * std::max(Scalar{1}, scale);
    auto predicate = [lhp_margin](std::complex<Scalar> lam) -> bool
    {
        return lam.real() < -lhp_margin;
    };

    auto rr = reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
    if (rr.placed < n)
        return std::unexpected(care_error::non_lhp_stabilisable);
    if (!T.allFinite() || !U.allFinite())
        return std::unexpected(care_error::non_finite_input);

    U.leftCols(n).array().colwise() *= D.array();

    care_result<Scalar, NX> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if (!P_err)
    {
        switch (P_err.error())
        {
            case riccati_extract_error::singular_u11:
                return std::unexpected(care_error::singular_u11);
            case riccati_extract_error::non_finite:
                return std::unexpected(care_error::non_finite_input);
            case riccati_extract_error::non_psd:
                return std::unexpected(care_error::non_psd_solution);
        }
        return std::unexpected(care_error::non_finite_input);
    }

    out.subspace_separation = rr.subspace_separation;
    out.reorder_complete    = rr.complete;
    return out;
}

}

#endif
