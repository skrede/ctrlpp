#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_H

/// @brief Discrete Algebraic Riccati Equation solver via real-Schur Bai-Demmel reorder.
///
/// Solves A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q = 0 for the stabilising P.
///
/// Builds the 2n x 2n symplectic matrix Z (Laub 1979 Eq. 7), computes its real Schur
/// decomposition Z = U T U^T, reorders T with a predicate `|lambda| < 1 - eps * scale`
/// via the Bai-Demmel 1993 swap kernel, and extracts P = U21 * U11^-1 from the resulting
/// invariant subspace basis. The conditioning of each adjacent-block swap is policed by
/// a compile-time policy (default `pivot_ratio_conditioning`).
///
/// @cite laub1979       -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
/// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include "ctrlpp/types.h"
#include "ctrlpp/control/dare_types.h"

#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <expected>
#include <type_traits>

namespace ctrlpp
{

namespace detail
{

/// @brief Build the symplectic matrix Z for DARE from (A, B, Q, R) per Laub 1979 Eq. 7.
///
/// Z = [[A + G A^{-T} Q,  -G A^{-T}],
///      [-A^{-T} Q,         A^{-T}  ]]   where G = B R^{-1} B^T
template <typename Scalar, std::size_t NX, std::size_t NU>
auto build_dare_symplectic(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                           const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                           const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                           const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> std::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, dare_error>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using MatNxN   = Eigen::Matrix<Scalar, n, n>;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;

    if (!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return std::unexpected(dare_error::non_finite_input);

    auto qr_At = A.transpose().colPivHouseholderQr();
    if (!qr_At.isInvertible())
        return std::unexpected(dare_error::non_finite_input);

    const MatNxN AinvT = qr_At.solve(MatNxN::Identity()).eval();
    const MatNxN G = (B * R.colPivHouseholderQr().solve(
                             Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval();

    Mat2Nx2N Z;
    Z.template block<n, n>(0, 0) = A + G * AinvT * Q;
    Z.template block<n, n>(0, n) = -G * AinvT;
    Z.template block<n, n>(n, 0) = -AinvT * Q;
    Z.template block<n, n>(n, n) = AinvT;

    if (!Z.allFinite())
        return std::unexpected(dare_error::non_finite_input);

    return Z;
}

}

/// @brief Discrete Algebraic Riccati Equation solver.
///
/// Returns `std::expected<dare_result<Scalar, NX>, dare_error>`. On success,
/// `result->P` is the stabilising solution; `result->subspace_separation` is the
/// min pivot ratio across accepted swaps (LAPACK SEP analogue); `result->reorder_complete`
/// is true iff every swap was accepted by the conditioning test.
template <typename Scalar, std::size_t NX, std::size_t NU,
          detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          Cond                                           /*tag*/ = {})
    -> std::expected<dare_result<Scalar, NX>, dare_error>
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    auto Z_result = detail::build_dare_symplectic<Scalar, NX, NU>(A, B, Q, R);
    if (!Z_result)
        return std::unexpected(Z_result.error());
    const Mat2N& Z = *Z_result;

    Eigen::RealSchur<Mat2N> schur(Z);
    if (schur.info() != Eigen::Success)
        return std::unexpected(dare_error::schur_failed);

    Mat2N T = schur.matrixT();
    Mat2N U = schur.matrixU();
    if (!T.allFinite() || !U.allFinite())
        return std::unexpected(dare_error::non_finite_input);

    const Scalar scale = T.cwiseAbs().maxCoeff();
    const Scalar eps   = std::numeric_limits<Scalar>::epsilon();
    const Scalar unit_margin = eps * std::max(Scalar{1}, scale);
    auto predicate = [unit_margin](std::complex<Scalar> lam) -> bool
    {
        return std::abs(lam) < Scalar{1} - unit_margin;
    };

    auto rr = detail::reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
    if (rr.placed < n)
        return std::unexpected(dare_error::non_stabilisable);
    if (!T.allFinite() || !U.allFinite())
        return std::unexpected(dare_error::non_finite_input);

    auto P_result = detail::extract_riccati_solution<Scalar, n2>(U);
    if (!P_result)
    {
        switch (P_result.error())
        {
            case detail::riccati_extract_error::singular_u11:
                return std::unexpected(dare_error::singular_u11);
            case detail::riccati_extract_error::non_finite:
                return std::unexpected(dare_error::non_finite_input);
            case detail::riccati_extract_error::non_psd:
                return std::unexpected(dare_error::non_psd_solution);
        }
        return std::unexpected(dare_error::non_finite_input);
    }

    dare_result<Scalar, NX> out;
    out.P                   = *P_result;
    out.subspace_separation = rr.subspace_separation;
    out.reorder_complete    = rr.complete;
    return out;
}

/// @brief DARE with cross-weight N: reduces to standard form via
/// Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T, then forwards.
template <typename Scalar, std::size_t NX, std::size_t NU,
          detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& N,
          Cond                                           tag = {})
    -> std::expected<dare_result<Scalar, NX>, dare_error>
{
    auto Rinv_Nt = R.colPivHouseholderQr()
                       .solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose()))
                       .eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return dare<Scalar, NX, NU, Cond>(Ap, B, Qp, R, tag);
}

}

#endif
