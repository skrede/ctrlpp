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
/// A small change between iterations is only a candidate for convergence: an
/// ill-conditioned iterate can stop changing relative to its norm before it is
/// a matrix sign. At a candidate, the stable-subspace projector
/// (I - H_k) / 2 is first checked for idempotence with a counted rounding
/// bound. If nonnormal rounding prevents that sufficient check from resolving,
/// the extracted solution must instead satisfy a counted Riccati residual
/// bound and place the closed-loop spectrum strictly in the open left
/// half-plane. These are properties of the fixed point and answer rather than
/// inferences from the iteration history. The stopping criterion is derived in
/// place; Higham supplies the scaled Newton iteration, not these acceptance
/// bounds. The rank threshold handed to Eigen is the relative, dimensionless
/// multiplier 2n times epsilon that its rank() compares against the largest
/// pivot. No bare numeric literal other than structural constants (the 1/2
/// from Eq. 5.16, the 2 from 2n = size(H), iteration cap 40 from Higham Table
/// 5.2's worst-observed count with determinantal scaling) appears in the hot
/// path.
///
/// @cite roberts1980 : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
/// @cite byers1987   : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
/// @cite higham2008  : Higham, "Functions of Matrices: Theory and Computation", 2008, Ch. 5 Eq. 5.16 / 5.34 / 5.35 / 5.41, Theorem 5.1(e)

#include "ctrlpp/expected.h"

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/QR>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
auto care_solution_satisfies_postconditions(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H,
    const Eigen::Matrix<Scalar, int(NX), int(NX)>& P) -> bool
{
    constexpr int n = int(NX);
    using MatN = Eigen::Matrix<Scalar, n, n>;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const MatN A = H.template block<n, n>(0, 0);
    const MatN S = -H.template block<n, n>(0, n);
    const MatN Q = -H.template block<n, n>(n, 0);
    const MatN At_P = (A.transpose() * P).eval();
    const MatN P_A = (P * A).eval();
    const MatN P_S_P = (P * S * P).eval();
    const MatN residual = (At_P + P_A - P_S_P + Q).eval();

    // The extracted P carries the arithmetic of its solve as well as the
    // residual evaluation. A conservative count uses 2n^3 operations for the
    // n-by-n Householder QR, 3n^3 for applying Q and solving n triangular
    // right-hand sides, and 3n^2 for symmetrization. Four n-by-n products then
    // use n^2 * (2n - 1) operations each, and three n-by-n sums combine the
    // residual terms. Every operation is counted whether or not it rounds.
    // The Frobenius scale is the largest term entering the final sum, not the
    // cancellation residual.
    constexpr int extraction_qr_rounding_ops =
        2 * n * n * n;
    constexpr int extraction_solve_rounding_ops =
        3 * n * n * n;
    constexpr int symmetrize_rounding_ops = 3 * n * n;
    constexpr int matrix_product_entry_rounding_ops = 2 * n - 1;
    constexpr int residual_matrix_products = 4;
    constexpr int residual_product_rounding_ops =
        residual_matrix_products * n * n
        * matrix_product_entry_rounding_ops;
    constexpr int residual_sums = 3;
    constexpr int residual_sum_rounding_ops =
        residual_sums * n * n;
    constexpr int residual_rounding_ops =
        extraction_qr_rounding_ops
        + extraction_solve_rounding_ops
        + symmetrize_rounding_ops
        + residual_product_rounding_ops
        + residual_sum_rounding_ops;
    const Scalar residual_scale = std::max(
        {At_P.norm(), P_A.norm(), P_S_P.norm(), Q.norm()});
    const Scalar residual_floor =
        Scalar{residual_rounding_ops} * eps * residual_scale;
    const bool residual_is_resolved =
        residual.norm() <= residual_floor;
    if (!residual_is_resolved)
        return false;

    const MatN closed_loop = (A - S * P).eval();
    Eigen::RealSchur<MatN> closed_loop_schur(closed_loop);
    if (closed_loop_schur.info() != Eigen::Success)
        return false;

    const MatN closed_loop_T = closed_loop_schur.matrixT();
    if (!closed_loop_T.allFinite())
        return false;

    // A backward-stable real Schur factor has an n-operation eigenvalue
    // uncertainty at the factor's largest-entry scale. Requiring every
    // diagonal real part below its negative bound certifies that the
    // closed-loop spectrum lies strictly in the open left half-plane.
    constexpr int closed_loop_eigenvalue_rounding_ops = n;
    const Scalar closed_loop_scale =
        closed_loop_T.cwiseAbs().maxCoeff();
    const Scalar closed_loop_margin =
        Scalar{closed_loop_eigenvalue_rounding_ops} * eps
        * closed_loop_scale;
    for (int index = 0; index < n; ++index)
    {
        if (!(closed_loop_T(index, index) < -closed_loop_margin))
            return false;
    }
    return true;
}

template <typename Scalar, std::size_t NX>
auto care_solve_via_sign_function(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in)
    -> ctrlpp::expected<care_result<Scalar, NX>, care_error>
{
    constexpr int n  = int(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    if (!H_in.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    Mat2N H = H_in;
    Mat2N H_prev;
    Mat2N P_LHP;
    bool projector_is_idempotent = false;

    const Scalar eps       = std::numeric_limits<Scalar>::epsilon();
    constexpr int max_iters = 40;
    // The change test is a gate for the fixed-point check, not an acceptance
    // decision. It avoids paying for an extra matrix product on iterations
    // that are still moving substantially.
    const Scalar change_candidate_dimension = Scalar{2} * Scalar{n2};
    const Scalar change_candidate_tolerance =
        std::sqrt(eps) * change_candidate_dimension;
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
            return ctrlpp::unexpected(care_error::sign_function_stagnated);

        const Mat2N H_inv = lu.solve(Mat2N::Identity()).eval();
        H = ((Scalar{1} / Scalar{2}) * (mu * H + (Scalar{1} / mu) * H_inv)).eval();

        if (!H.allFinite())
            return ctrlpp::unexpected(care_error::non_finite_input);

        const Scalar delta_norm = (H - H_prev).norm();
        const Scalar scale_H    = H.norm();
        const Scalar change_candidate_floor =
            change_candidate_tolerance * scale_H;
        const bool change_is_small = delta_norm <= change_candidate_floor;
        if (change_is_small)
        {
            P_LHP = ((Scalar{1} / Scalar{2})
                     * (Mat2N::Identity() - H)).eval();
            const Mat2N projector_squared = (P_LHP * P_LHP).eval();
            const Scalar projector_defect =
                (projector_squared - P_LHP).norm();

            // Each entry of P^2 uses n2 products and n2 - 1 additions;
            // subtracting P adds one more rounding. That is 2 * n2 rounded
            // operations. The Frobenius operand norm below already aggregates
            // the entries, so the operation count is not multiplied by the
            // number of entries a second time. The operand scale is the larger
            // norm of the two matrices entering the final subtraction, not the
            // cancellation residual produced by it.
            constexpr int projector_idempotence_rounding_ops = 2 * n2;
            const Scalar projector_idempotence_scale =
                std::max(projector_squared.norm(), P_LHP.norm());
            const Scalar projector_idempotence_floor =
                Scalar{projector_idempotence_rounding_ops} * eps
                * projector_idempotence_scale;
            projector_is_idempotent =
                projector_defect <= projector_idempotence_floor;
            break;
        }

        // A change that grows after the warm-up window is non-contraction. A
        // small change cannot bypass the fixed-point check above, so this guard
        // remains separate from the resolved-sign decision.
        if (k > 3 && delta_norm > (Scalar{1} / Scalar{2}) * last_delta_norm)
            return ctrlpp::unexpected(care_error::sign_function_stagnated);
        last_delta_norm = delta_norm;

        if (k + 1 == max_iters)
            return ctrlpp::unexpected(care_error::sign_function_stagnated);
    }

    Eigen::ColPivHouseholderQR<Mat2N> qr(P_LHP);
    // Eigen's ColPivHouseholderQR::rank() compares each pivot against
    // threshold() times the largest pivot, so setThreshold takes a relative,
    // dimensionless multiplier. The backward-stable rank tolerance for a QR is
    // the matrix size times unit roundoff (2n times epsilon); multiplying by an
    // operand norm would apply the scale twice and inflate the cutoff on
    // exactly the ill-conditioned stable subspaces this path must resolve.
    qr.setThreshold(Scalar{n2} * eps);
    if (qr.rank() < n)
        return ctrlpp::unexpected(care_error::non_lhp_stabilisable);

    const Mat2N Q_full = qr.householderQ();
    const Mat2N U      = Q_full;

    care_result<Scalar, NX> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if (!P_err)
    {
        switch (P_err.error())
        {
            case riccati_extract_error::singular_u11:
                return ctrlpp::unexpected(care_error::singular_u11);
            case riccati_extract_error::non_finite:
                return ctrlpp::unexpected(care_error::non_finite_input);
            case riccati_extract_error::non_psd:
                return ctrlpp::unexpected(care_error::non_psd_solution);
        }
        return ctrlpp::unexpected(care_error::non_finite_input);
    }

    if (!projector_is_idempotent
        && !care_solution_satisfies_postconditions<Scalar, NX>(
            H_in, out.P))
        return ctrlpp::unexpected(
            care_error::sign_function_stagnated);

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
