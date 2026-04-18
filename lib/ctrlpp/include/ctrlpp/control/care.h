#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_H

/// @brief Continuous-time Algebraic Riccati Equation solver via real-Schur Bai-Demmel reorder.
///
/// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0 for the stabilising P.
///
/// Builds the 2n x 2n Hamiltonian H = [[A, -B R^{-1} B^T], [-Q, -A^T]] (Laub 1979),
/// computes its real Schur decomposition H = U T U^T, reorders T with a predicate
/// `Re(lambda) < -eps * scale` (open left half-plane) via the Bai-Demmel 1993 swap
/// kernel, and extracts P = U21 * U11^-1 from the resulting invariant subspace basis.
///
/// The continuous-time LQR gain is K = R^{-1} B^T P.
///
/// Shares the real-Schur reorder primitive and the Riccati P-extraction helper with
/// `ctrlpp::dare` via `ctrlpp::detail/`. Does not depend on `dare.h`.
///
/// @cite laub1979       -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
/// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include "ctrlpp/types.h"
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
#include <type_traits>

namespace ctrlpp
{

namespace detail
{

/// @brief Build the 2n x 2n Hamiltonian matrix H for CARE from (A, B, Q, R).
///
/// H = [[A,    -B R^{-1} B^T],
///      [-Q,   -A^T         ]]
template <typename Scalar, std::size_t NX, std::size_t NU>
auto build_care_hamiltonian(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                            const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                            const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                            const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> std::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, care_error>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using MatNxN   = Eigen::Matrix<Scalar, n, n>;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;

    if (!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return std::unexpected(care_error::non_finite_input);

    const MatNxN S = (B * R.colPivHouseholderQr().solve(
                             Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval();

    Mat2Nx2N H;
    H.template block<n, n>(0, 0) = A;
    H.template block<n, n>(0, n) = -S;
    H.template block<n, n>(n, 0) = -Q;
    H.template block<n, n>(n, n) = -A.transpose();

    if (!H.allFinite())
        return std::unexpected(care_error::non_finite_input);

    return H;
}

}

/// @brief Continuous-time Algebraic Riccati Equation solver.
///
/// Returns `std::expected<care_result<Scalar, NX>, care_error>`. On success,
/// `result->P` is the stabilising solution; `result->subspace_separation` is the
/// min pivot ratio across accepted swaps; `result->reorder_complete` is true iff
/// every swap was accepted.
template <typename Scalar, std::size_t NX, std::size_t NU,
          detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          Cond                                           /*tag*/ = {})
    -> std::expected<care_result<Scalar, NX>, care_error>
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    auto H_result = detail::build_care_hamiltonian<Scalar, NX, NU>(A, B, Q, R);
    if (!H_result)
        return std::unexpected(H_result.error());
    const Mat2N& H = *H_result;

    Eigen::RealSchur<Mat2N> schur(H);
    if (schur.info() != Eigen::Success)
        return std::unexpected(care_error::schur_failed);

    Mat2N T = schur.matrixT();
    Mat2N U = schur.matrixU();
    if (!T.allFinite() || !U.allFinite())
        return std::unexpected(care_error::non_finite_input);

    const Scalar scale = T.cwiseAbs().maxCoeff();
    const Scalar eps   = std::numeric_limits<Scalar>::epsilon();
    const Scalar lhp_margin = eps * std::max(Scalar{1}, scale);
    auto predicate = [lhp_margin](std::complex<Scalar> lam) -> bool
    {
        return lam.real() < -lhp_margin;
    };

    auto rr = detail::reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
    if (rr.placed < n)
        return std::unexpected(care_error::non_lhp_stabilisable);
    if (!T.allFinite() || !U.allFinite())
        return std::unexpected(care_error::non_finite_input);

    auto P_result = detail::extract_riccati_solution<Scalar, n2>(U);
    if (!P_result)
    {
        switch (P_result.error())
        {
            case detail::riccati_extract_error::singular_u11:
                return std::unexpected(care_error::singular_u11);
            case detail::riccati_extract_error::non_finite:
                return std::unexpected(care_error::non_finite_input);
            case detail::riccati_extract_error::non_psd:
                return std::unexpected(care_error::non_psd_solution);
        }
        return std::unexpected(care_error::non_finite_input);
    }

    care_result<Scalar, NX> out;
    out.P                   = *P_result;
    out.subspace_separation = rr.subspace_separation;
    out.reorder_complete    = rr.complete;
    return out;
}

/// @brief CARE with cross-weight N: reduces to standard form via
/// Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T, then forwards.
template <typename Scalar, std::size_t NX, std::size_t NU,
          detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& N,
          Cond                                           tag = {})
    -> std::expected<care_result<Scalar, NX>, care_error>
{
    auto Rinv_Nt = R.colPivHouseholderQr()
                       .solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose()))
                       .eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return care<Scalar, NX, NU, Cond>(Ap, B, Qp, R, tag);
}

}

#endif
