#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_SIGN_FUNCTION_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_SIGN_FUNCTION_H

/// @brief Continuous-time Riccati solve via the matrix sign function.
///
/// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0 for the stabilising P by
/// iterating the Roberts 1980 Newton form
///     H_{k+1} = 0.5 * (mu_k * H_k + mu_k^{-1} * H_k^{-1})
/// on the Hamiltonian H until it converges to sign(H). The determinantal
/// scaling factor mu_k = |det(H_k)|^{-1/(2n)} is computed in overflow-safe
/// form from the LU factor diagonal per Higham 2008 Sec. 5.5:
///     mu_k = exp(-(1/(2n)) * sum_i log|u_ii|)
/// where u_ii are the diagonal entries of the Eigen::PartialPivLU matrixLU
/// already produced for the Newton-step inverse. The stable invariant
/// subspace is the range of the projector (I - sign(H))/2 (Higham Theorem
/// 5.1(e)); a rank-revealing Eigen::ColPivHouseholderQR on that projector
/// yields the 2n x n LHP subspace basis in the leading n columns of Q. The
/// Riccati solution is then P = U21 * U11^{-1} via the shared
/// ctrlpp::detail::extract_riccati_solution_into primitive, matching the
/// convention of the Schur-based path.
///
/// The Newton stopping test scales the square root of
/// std::numeric_limits<Scalar>::epsilon() by a dimensionless factor from the
/// size 2n: the iteration converges quadratically, so a relative change of
/// sqrt(epsilon) leaves the current iterate accurate to epsilon, and a
/// linear-in-epsilon gate would be unreachable on an ill-conditioned
/// Hamiltonian before its rounding floor. The rank threshold handed to Eigen
/// is the relative, dimensionless multiplier 2n times epsilon that its rank()
/// compares against the largest pivot. No bare numeric literal other than
/// structural constants (the 1/2 from Eq. 5.16, the 2 from 2n = size(H),
/// iteration cap 40 from Higham Table 5.2's worst-observed count with
/// determinantal scaling) appears in the hot path.
///
/// @cite roberts1980 : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
/// @cite byers1987   : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
/// @cite higham2008  : Higham, "Functions of Matrices: Theory and Computation", 2008, Ch. 5 Eq. 5.16 / 5.34 / 5.35 / 5.41, Theorem 5.1(e)

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>
#include <expected>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
auto care_solve_via_sign_function(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in)
    -> std::expected<care_result<Scalar, NX>, care_error>
{
    constexpr int n  = int(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    if (!H_in.allFinite())
        return std::unexpected(care_error::non_finite_input);

    Mat2N H = H_in;
    Mat2N H_prev;

    const Scalar eps       = std::numeric_limits<Scalar>::epsilon();
    constexpr int max_iters = 40;
    // The Newton sign iteration converges quadratically near the fixed point,
    // so successive changes satisfy delta_{k+1} ~ delta_k^2: once the relative
    // change reaches sqrt(epsilon), the current iterate already carries the
    // full epsilon-level accuracy of sign(H). The stopping test therefore
    // compares the relative change against sqrt(epsilon) times a dimensionless
    // factor from the matrix size 2n, which stays reachable on ill-conditioned
    // Hamiltonians where a linear-in-epsilon gate would sit below the rounding
    // floor (Higham 2008 Sec. 5.5, scaled-Newton stopping test).
    const Scalar conv_dim  = Scalar{2} * Scalar{n2};
    const Scalar conv_tol  = std::sqrt(eps) * conv_dim;
    Scalar last_delta_norm = std::numeric_limits<Scalar>::infinity();

    for (int k = 0; k < max_iters; ++k)
    {
        H_prev = H;

        Eigen::PartialPivLU<Mat2N> lu(H);
        Scalar log_det_abs = Scalar{0};
        for (int i = 0; i < n2; ++i)
            log_det_abs += std::log(std::abs(lu.matrixLU()(i, i)));
        const Scalar mu = std::exp(-log_det_abs / Scalar{n2});

        if (!std::isfinite(mu))
            return std::unexpected(care_error::sign_function_stagnated);

        const Mat2N H_inv = lu.solve(Mat2N::Identity()).eval();
        H = ((Scalar{1} / Scalar{2}) * (mu * H + (Scalar{1} / mu) * H_inv)).eval();

        if (!H.allFinite())
            return std::unexpected(care_error::non_finite_input);

        const Scalar delta_norm = (H - H_prev).norm();
        const Scalar scale_H    = H.norm();
        if (delta_norm <= conv_tol * scale_H)
            break;

        // Genuine non-monotone growth: the relative change is increasing while
        // still above the sqrt(epsilon) convergence plateau, which indicates
        // divergence rather than settling. The break above already accepts an
        // iterate that has reached the quadratic-convergence plateau, so this
        // guard cannot fire on a converged solve merely sitting at the floor.
        if (k > 3 && delta_norm > (Scalar{1} / Scalar{2}) * last_delta_norm)
            return std::unexpected(care_error::sign_function_stagnated);
        last_delta_norm = delta_norm;

        if (k + 1 == max_iters)
            return std::unexpected(care_error::sign_function_stagnated);
    }

    const Mat2N P_LHP = ((Scalar{1} / Scalar{2}) * (Mat2N::Identity() - H)).eval();

    Eigen::ColPivHouseholderQR<Mat2N> qr(P_LHP);
    // Eigen's ColPivHouseholderQR::rank() compares each pivot against
    // threshold() times the largest pivot, so setThreshold takes a relative,
    // dimensionless multiplier. The backward-stable rank tolerance for a QR is
    // the matrix size times unit roundoff (2n times epsilon); multiplying by an
    // operand norm would apply the scale twice and inflate the cutoff on
    // exactly the ill-conditioned stable subspaces this path must resolve.
    qr.setThreshold(Scalar{n2} * eps);
    if (qr.rank() < n)
        return std::unexpected(care_error::non_lhp_stabilisable);

    const Mat2N Q_full = qr.householderQ();
    const Mat2N U      = Q_full;

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

    // The sign-function path has no swap phase and hence no rank-revealing QR pivot ratio
    // to report; writing quiet_NaN signals "unavailable for this method" per the
    // care_result::subspace_separation contract, in contrast to the Schur path which
    // reports a finite positive pivot ratio.
    out.subspace_separation = std::numeric_limits<Scalar>::quiet_NaN();
    out.reorder_complete    = true;
    return out;
}

}

#endif
