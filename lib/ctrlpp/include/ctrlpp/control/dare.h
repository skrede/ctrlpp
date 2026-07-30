#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_H

/// @brief Discrete Algebraic Riccati Equation solver via real-Schur Bai-Demmel reorder.
///
/// Solves A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q = 0 for the stabilizing P.
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
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/dare_types.h"

#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <algorithm>

namespace ctrlpp {

namespace detail {

template<typename Scalar, std::size_t NX>
struct dare_symplectic_operands
{
    Eigen::Matrix<Scalar, int(NX), int(NX)> A_inverse_transpose;
    Eigen::Matrix<Scalar, int(NX), int(NX)> G;
};

template<typename Scalar, std::size_t NX, std::size_t NU>
auto factor_dare_g(const Eigen::Matrix<Scalar, int(NX), int(NU)> &B, const Eigen::Matrix<Scalar, int(NU), int(NU)> &R)
        -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NX), int(NX)>, dare_error>
{
    // G requires a genuine inverse of R. A rank-deficient factor would make
    // solve() return a least-squares result and silently pose a different
    // control problem, so apply a dimension-scaled reciprocal-pivot test first.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_R.isInvertible())
        return ctrlpp::unexpected(dare_error::singular_r);

    return Eigen::Matrix<Scalar, int(NX), int(NX)>{(B * qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval()};
}

/// @brief Factor the scale-invariant operands used by the DARE symplectic matrix.
template<typename Scalar, std::size_t NX, std::size_t NU>
auto factor_dare_symplectic_operands(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B,
                                     const Eigen::Matrix<Scalar, int(NU), int(NU)> &R) -> ctrlpp::expected<dare_symplectic_operands<Scalar, NX>, dare_error>
{
    constexpr int n = static_cast<int>(NX);
    using MatNxN    = Eigen::Matrix<Scalar, n, n>;

    // A^{-T} exists only when A has full rank. The state dimension times unit
    // roundoff supplies the scale-relative reciprocal-pivot threshold.
    auto qr_At = A.transpose().colPivHouseholderQr();
    qr_At.setThreshold(Scalar{static_cast<int>(NX)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_At.isInvertible())
        return ctrlpp::unexpected(dare_error::singular_a);

    auto G_result = factor_dare_g<Scalar, NX, NU>(B, R);
    if(!G_result)
        return ctrlpp::unexpected(G_result.error());

    dare_symplectic_operands<Scalar, NX> operands;
    operands.A_inverse_transpose = qr_At.solve(MatNxN::Identity()).eval();
    operands.G                   = *G_result;
    return operands;
}

/// @brief Choose the common weight divisor from the largest weight entry.
///
/// Dividing Q and R by the same positive value preserves the Riccati gain.
/// Selecting max(||Q||_max, ||R||_max) gives the equivalent posed problem a
/// canonical scale and keeps at least one weight entry near unity.
template<typename Scalar, std::size_t NX, std::size_t NU>
auto dare_weight_scale(const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q, const Eigen::Matrix<Scalar, int(NU), int(NU)> &R) -> Scalar
{
    const Scalar q_magnitude      = Q.cwiseAbs().maxCoeff();
    const Scalar r_magnitude      = R.cwiseAbs().maxCoeff();
    const Scalar weight_magnitude = std::max(q_magnitude, r_magnitude);
    return weight_magnitude > Scalar{0} ? weight_magnitude : Scalar{1};
}

/// @brief Build the equilibrated symplectic matrix Z per Laub 1979 Eq. 7.
///
/// Z = [[A + G A^{-T} Q,  -G A^{-T}],
///      [-A^{-T} Q,         A^{-T}  ]]   where G = B R^{-1} B^T
template<typename Scalar, std::size_t NX>
auto build_dare_symplectic(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q, const dare_symplectic_operands<Scalar, NX> &operands)
        -> ctrlpp::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, dare_error>
{
    constexpr int n  = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using Mat2Nx2N   = Eigen::Matrix<Scalar, n2, n2>;

    Mat2Nx2N Z;
    Z.template block<n, n>(0, 0) = A + operands.G * operands.A_inverse_transpose * Q;
    Z.template block<n, n>(0, n) = -operands.G * operands.A_inverse_transpose;
    Z.template block<n, n>(n, 0) = -operands.A_inverse_transpose * Q;
    Z.template block<n, n>(n, n) = operands.A_inverse_transpose;

    if(!Z.allFinite())
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    return Z;
}

// Six state-dimension contractions enter A'PA, B'PA and A'PBK; two
// input-dimension contractions enter B'PB and the final gain product.
// One sum forms R + B'PB, the rank-revealing gain solve contributes two
// input-dimension operations, and three sums assemble the four residual
// terms. Each length-d contraction has d multiplies and d-1 additions.
template<std::size_t NX, std::size_t NU>
constexpr int dare_residual_ops = 6 * (2 * static_cast<int>(NX) - 1) + 2 * (2 * static_cast<int>(NU) - 1) + 1 + 2 * static_cast<int>(NU) + 3;

template<typename Scalar, std::size_t NX, std::size_t NU>
auto compute_dare_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B, const Eigen::Matrix<Scalar, int(NU), int(NU)> &R,
                       const Eigen::Matrix<Scalar, int(NX), int(NX)> &P, Eigen::Matrix<Scalar, int(NU), int(NX)> &gain) -> bool
{
    const auto BtP          = (B.transpose() * P).eval();
    const auto gain_operand = (R + BtP * B).eval();
    auto gain_qr            = gain_operand.colPivHouseholderQr();
    gain_qr.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!gain_qr.isInvertible())
        return false;

    gain = gain_qr.solve(BtP * A).eval();
    return gain.allFinite();
}

/// @brief Verify that P solves the posed DARE and produces a stabilizing gain.
template<typename Scalar, std::size_t NX, std::size_t NU>
auto verify_dare_solution(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B, const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q,
                          const Eigen::Matrix<Scalar, int(NU), int(NU)> &R, const Eigen::Matrix<Scalar, int(NX), int(NX)> &P,
                          Eigen::Matrix<Scalar, int(NU), int(NX)> *verified_gain = nullptr) -> bool
{
    constexpr int n = static_cast<int>(NX);
    using MatNxN    = Eigen::Matrix<Scalar, n, n>;

    Eigen::Matrix<Scalar, int(NU), int(NX)> gain;
    if(!compute_dare_gain<Scalar, NX, NU>(A, B, R, P, gain))
        return false;

    const MatNxN AtPA     = (A.transpose() * P * A).eval();
    const MatNxN AtPBK    = (A.transpose() * P * B * gain).eval();
    const MatNxN residual = (AtPA - P - AtPBK + Q).eval();
    if(!residual.allFinite())
        return false;

    const Scalar residual_scale     = std::max({AtPA.norm(), P.norm(), AtPBK.norm(), Q.norm()});
    const Scalar residual_magnitude = residual.norm();
    if(!std::isfinite(residual_scale) || !std::isfinite(residual_magnitude))
        return false;
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    // A Schur-extracted invariant subspace is a forward solution, not a direct
    // evaluation of the residual expression. Requiring its forward residual to
    // sit at the expression's pure rounding floor over-refuses correctly solved
    // conditioned problems. The square root of the enumerated rounding budget
    // is the precision boundary: it retains at least half the scalar type's
    // significand while still declining the measured wrong-success population.
    const Scalar residual_margin = std::sqrt(Scalar{dare_residual_ops<NX, NU>} * eps) * residual_scale;
    if(!(residual_magnitude <= residual_margin))
        return false;

    const MatNxN closed_loop = (A - B * gain).eval();
    if(!closed_loop.allFinite())
        return false;
    Eigen::EigenSolver<MatNxN> eigensolver(closed_loop, false);
    if(eigensolver.info() != Eigen::Success || !eigensolver.eigenvalues().allFinite())
        return false;

    const Scalar closed_loop_scale = std::max(Scalar{1}, closed_loop.cwiseAbs().maxCoeff());
    // Each closed-loop eigenvalue carries the standard n * eps * ||A-BK||
    // backward-error margin. A pole inside the disk by less than this amount
    // cannot support the solver's stabilizing claim at the represented scale.
    const Scalar unit_margin = Scalar{static_cast<int>(NX)} * eps * closed_loop_scale;
    for(int index = 0; index < n; ++index)
    {
        if(!(std::abs(eigensolver.eigenvalues()(index)) < Scalar{1} - unit_margin))
            return false;
    }
    if(verified_gain != nullptr)
        *verified_gain = gain;
    return true;
}

/// @brief Solve DARE from a pre-built symplectic Z: real-Schur, Bai-Demmel reorder inside
/// the unit disk, Riccati extract.
template<typename Scalar, std::size_t NX, conditioning_policy Cond = pivot_ratio_conditioning>
auto dare_solve_from_symplectic(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)> &Z, Cond /*tag*/ = {}) -> ctrlpp::expected<dare_result<Scalar, NX>, dare_error>
{
    constexpr int n  = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using Mat2N      = Eigen::Matrix<Scalar, n2, n2>;

    Eigen::RealSchur<Mat2N> schur(Z);
    if(schur.info() != Eigen::Success)
        return ctrlpp::unexpected(dare_error::schur_failed);

    Mat2N T = schur.matrixT();
    Mat2N U = schur.matrixU();
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    const Scalar scale = T.cwiseAbs().maxCoeff();
    const Scalar eps   = std::numeric_limits<Scalar>::epsilon();
    // Eigenvalues of a backward-stable real Schur factor carry a perturbation
    // on the order of the matrix size times unit roundoff times the factor
    // norm, so the unit-disk predicate margin is that backward error: 2n times
    // epsilon times the largest magnitude of T.
    const Scalar unit_margin = Scalar{n2} * eps * scale;
    auto predicate           = [unit_margin](std::complex<Scalar> lam) -> bool { return std::abs(lam) < Scalar{1} - unit_margin; };

    auto rr = reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
    if(rr.placed < n)
        return ctrlpp::unexpected(dare_error::non_stabilizable);
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    dare_result<Scalar, NX> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if(!P_err)
    {
        switch(P_err.error())
        {
            case riccati_extract_error::singular_u11:
                return ctrlpp::unexpected(dare_error::singular_u11);
            case riccati_extract_error::non_finite:
                return ctrlpp::unexpected(dare_error::non_finite_input);
            case riccati_extract_error::non_psd:
                return ctrlpp::unexpected(dare_error::non_psd_solution);
        }
        return ctrlpp::unexpected(dare_error::non_finite_input);
    }

    out.subspace_separation = rr.subspace_separation;
    out.reorder_complete    = rr.complete;
    return out;
}

}

/// @brief Discrete Algebraic Riccati Equation solver.
///
/// Returns `ctrlpp::expected<dare_result<Scalar, NX>, dare_error>`. On success,
/// `result->P` is the stabilizing solution; `result->subspace_separation` is the
/// min pivot ratio across accepted swaps (LAPACK SEP analogue); `result->reorder_complete`
/// is true iff every swap was accepted by the conditioning test.
template<ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B, const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)> &R, Cond /*tag*/ = {}) -> ctrlpp::expected<dare_result<Scalar, NX>, dare_error>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    auto operands_result = detail::factor_dare_symplectic_operands<Scalar, NX, NU>(A, B, R);
    if(!operands_result)
        return ctrlpp::unexpected(operands_result.error());

    const Scalar weight_scale = detail::dare_weight_scale<Scalar, NX, NU>(Q, R);
    if(!(weight_scale > Scalar{0}) || !std::isfinite(weight_scale))
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    auto Q_scaled        = Q;
    auto R_scaled        = R;
    auto scaled_operands = *operands_result;
    if(weight_scale != Scalar{1})
    {
        Q_scaled /= weight_scale;
        R_scaled /= weight_scale;
        if(!Q_scaled.allFinite() || !R_scaled.allFinite())
            return ctrlpp::unexpected(dare_error::arithmetic_limit);

        auto scaled_G = detail::factor_dare_g<Scalar, NX, NU>(B, R_scaled);
        if(!scaled_G)
            return ctrlpp::unexpected(dare_error::arithmetic_limit);
        scaled_operands.G = *scaled_G;
    }

    auto Z_result = detail::build_dare_symplectic<Scalar, NX>(A, Q_scaled, scaled_operands);
    if(!Z_result)
        return ctrlpp::unexpected(Z_result.error());

    auto result = detail::dare_solve_from_symplectic<Scalar, NX, Cond>(*Z_result);
    if(!result)
    {
        if(result.error() == dare_error::non_finite_input)
            return ctrlpp::unexpected(dare_error::arithmetic_limit);
        return ctrlpp::unexpected(result.error());
    }

    Eigen::Matrix<Scalar, int(NU), int(NX)> scaled_gain;
    if(!Q_scaled.allFinite() || !R_scaled.allFinite() || !detail::verify_dare_solution<Scalar, NX, NU>(A, B, Q_scaled, R_scaled, result->P, &scaled_gain))
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    if(weight_scale == Scalar{1})
        return result;

    result->P *= weight_scale;
    Eigen::Matrix<Scalar, int(NU), int(NX)> returned_gain;
    if(!result->P.allFinite() || !detail::verify_dare_solution<Scalar, NX, NU>(A, B, Q, R, result->P, &returned_gain))
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    const Scalar gain_scale  = std::max(scaled_gain.norm(), returned_gain.norm());
    const Scalar gain_margin = std::sqrt(Scalar{detail::dare_residual_ops<NX, NU>} * std::numeric_limits<Scalar>::epsilon()) * gain_scale;
    if(!((scaled_gain - returned_gain).norm() <= gain_margin))
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    return result;
}

/// @brief DARE with cross-weight N: reduces to standard form via
/// Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T, then forwards.
template<typename Scalar, std::size_t NX, std::size_t NU, detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B, const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)> &R, const Eigen::Matrix<Scalar, int(NX), int(NU)> &N, Cond tag = {}) -> ctrlpp::expected<dare_result<Scalar, NX>, dare_error>
{
    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite() || !N.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    // The reduction to standard form inverts R before the symplectic build ever sees
    // it, so the same test is applied to the factorization this overload already forms
    // rather than deferring to the build's. Without it a singular R would reach the
    // inner solve as a non-finite Q' and A' and be reported as a non-finite input.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_R.isInvertible())
        return ctrlpp::unexpected(dare_error::singular_r);

    auto Rinv_Nt = qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose())).eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();
    if(!Qp.allFinite() || !Ap.allFinite())
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    return dare<Scalar, NX, NU, Cond>(Ap, B, Qp, R, tag);
}

}

#endif
