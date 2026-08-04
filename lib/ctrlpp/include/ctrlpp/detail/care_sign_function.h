#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_SIGN_FUNCTION_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_SIGN_FUNCTION_H

/// @brief Continuous-time Riccati solve via the matrix sign function.
///
/// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0 for the stabilizing P by
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
/// a matrix sign. Acceptance is therefore decided on the RETURNED solution and
/// on nothing else: the extracted P must satisfy a counted Riccati residual
/// bound and place the closed-loop spectrum strictly in the open left
/// half-plane. These are properties of the answer rather than inferences from
/// the iteration history.
///
/// That rule is not this path's private property. It lives in
/// ctrlpp/detail/care_postconditions.h and is the SAME rule the real-Schur and
/// balanced-Schur paths reach, because the public result contract promises a
/// stabilizing matrix without qualifying by method. What remains specific to
/// this path is only which enumerator reports a failure of it: an iteration
/// that produced an answer this path could not verify is reported as
/// stagnation, alongside the other three ways the iteration itself can fail.
///
/// The projector's idempotence is deliberately NOT part of that decision, and
/// this is a correction rather than a simplification. Idempotence of
/// (I - H_k) / 2 is evidence that the iteration reached a fixed point; it is
/// not evidence about the matrix extracted from that fixed point. Measurement
/// separates the two: across the whole magnitude band in which the extracted
/// solution comes back as exactly zero, the projector's idempotence defect is
/// exactly zero against its rounding floor at every point, while the residual
/// postcondition on the returned matrix is false at every point. A verdict that
/// cannot distinguish the two cannot be allowed to certify one of them, so the
/// postconditions run unconditionally and the redundant projector product is
/// gone from the path.
///
/// The stopping criterion is derived in place; Higham supplies the scaled
/// Newton iteration, not these acceptance bounds. The rank threshold handed to
/// Eigen is the relative, dimensionless multiplier 2n times epsilon that its
/// rank() compares against the largest pivot. Every acceptance comparison is
/// evaluated on ctrlpp::detail::resolved_magnitude operands so that neither end
/// of the arithmetic range can turn a comparison into an unconditional verdict.
/// No bare numeric literal other than structural constants (the 1/2 from
/// Eq. 5.16, the 2 from 2n = size(H), iteration cap 40 from Higham Table 5.2's
/// worst-observed count with determinantal scaling) appears in the hot path.
///
/// @cite roberts1980 : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
/// @cite byers1987   : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
/// @cite higham2008  : Higham, "Functions of Matrices: Theory and Computation", 2008, Ch. 5 Eq. 5.16 / 5.34 / 5.35 / 5.41, Theorem 5.1(e)

#include "ctrlpp/expected.h"

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/care_postconditions.h"

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp::detail
{

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

    const Scalar eps       = std::numeric_limits<Scalar>::epsilon();
    constexpr int max_iters = 40;
    // The change test decides when to STOP iterating; it decides nothing about
    // the answer. Leaving the loop leads to the extraction and to the
    // verification of what was extracted, and it is that verification, not this
    // test, that permits a success.
    const Scalar change_candidate_dimension = Scalar{2} * Scalar{n2};
    const Scalar change_candidate_tolerance =
        std::sqrt(eps) * change_candidate_dimension;
    resolved_magnitude<Scalar> last_change{Scalar{0}, false};

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

        // On the first step the previous iterate is the caller's own
        // Hamiltonian, whose sum of squares can be past the top of the range
        // while every entry of it is finite. Resolving both operands is what
        // keeps this a comparison rather than a formality, and an operand that
        // does not resolve leaves through the typed channel instead of being
        // carried forward.
        const resolved_magnitude<Scalar> change =
            resolve_magnitude(H - H_prev);
        const resolved_magnitude<Scalar> scale_H = resolve_magnitude(H);
        if (!change.resolved || !scale_H.resolved)
            return ctrlpp::unexpected(care_error::sign_function_stagnated);

        const resolved_magnitude<Scalar> change_candidate_floor =
            scaled_magnitude(scale_H, change_candidate_tolerance);
        if (magnitude_within(change, change_candidate_floor))
        {
            P_LHP = ((Scalar{1} / Scalar{2})
                     * (Mat2N::Identity() - H)).eval();
            break;
        }

        // A change that grows after the warm-up window is non-contraction, and
        // an iteration that is not contracting will not reach a sign matrix.
        // This is a statement about the ITERATION, which is why it is separate
        // from the verification of the answer: it declines early rather than
        // spending the remaining budget to decline later.
        if (k > 3
            && !magnitude_within(
                change,
                scaled_magnitude(last_change,
                                 Scalar{1} / Scalar{2})))
            return ctrlpp::unexpected(care_error::sign_function_stagnated);
        last_change = change;

        if (k + 1 == max_iters)
            return ctrlpp::unexpected(care_error::sign_function_stagnated);
    }

    // The rank-revealing QR below is equivariant under a positive scalar
    // multiple of its operand: the orthogonal factor is unchanged, the
    // triangular factor and every column norm scale together so the pivot order
    // is unchanged, and a power-of-two factor is exact in binary floating point.
    // Rescaling the projector therefore changes no mathematics.
    //
    // It changes what the factorization can see. Each Householder step compares
    // the SQUARED norm of its column tail against an absolute floor -- the
    // scalar type's smallest normal value -- and discards the reflection when
    // the tail falls to or below it. The tail carries the projector's
    // stable/unstable coupling, whose size is the size of the solution, so a
    // solution smaller than the square root of the smallest normal value had its
    // reflection discarded and came back as exactly zero. Scaling the operand by
    // s scales every tail's squared norm by s^2 and lifts it clear.
    //
    // The factor is derived from the operand rather than chosen. A column tail
    // holds at most n2 entries, so bounding the operand's largest entry by
    // sqrt(max / n2) bounds every tail's sum of squares by the largest finite
    // value. std::frexp writes the ratio as m * 2^e with m in [0.5, 1), so the
    // largest power of two at or below it is 2^(e - 1), and that choice meets
    // the bound EXACTLY when the ratio is itself a power of two -- a sum of
    // squares landing on the largest finite value can round past it. One binary
    // order below, 2^(e - 2), bounds every tail by a quarter of the largest
    // finite value instead, so the bound is met strictly. The order costs a
    // factor of four in tail squared norm at the bottom, which is not the
    // binding constraint there: what is left, s^2, is still the whole headroom
    // the scalar type has, and the reach is set by the Newton step below rather
    // than by this factorization.
    //
    // The rescale is applied unconditionally and the postconditions decide what
    // it recovered. Its reach is not the QR's any more: past one over the square
    // root of the smallest normal value the Newton step's own inverse carries a
    // reciprocal squared magnitude into the subnormal range, and past one over
    // the square root of the smallest subnormal value that entry is exactly zero
    // and the accepted iterate carries exactly half the true coupling. No
    // rescale of the projector restores an answer the iterate no longer holds,
    // so the boundary of the recovered region is the residual postcondition
    // itself rather than a constant.
    Scalar projector_rescale = Scalar{1};
    const resolved_magnitude<Scalar> projector_largest_entry =
        resolve_largest_entry(P_LHP);
    if (projector_largest_entry.resolved
        && projector_largest_entry.value > Scalar{0})
    {
        const Scalar tail_ceiling =
            std::sqrt(std::numeric_limits<Scalar>::max() / Scalar{n2});
        const Scalar headroom_ratio =
            tail_ceiling / projector_largest_entry.value;
        if (std::isfinite(headroom_ratio) && headroom_ratio > Scalar{0})
        {
            int headroom_exponent = 0;
            std::frexp(headroom_ratio, &headroom_exponent);
            const Scalar candidate =
                std::ldexp(Scalar{1}, headroom_exponent - 2);
            if (std::isfinite(candidate) && candidate > Scalar{0})
                projector_rescale = candidate;
        }
    }

    // This is the factorization whose orthogonal factor becomes the basis the
    // extraction reads, and it is 2n by 2n. The shared acceptance rule counts
    // it at that size through care_extraction_basis_dimension rather than at
    // the n by n block the extraction writes into, so the residual bound this
    // path is held to covers the arithmetic this path actually performed.
    Eigen::ColPivHouseholderQR<Mat2N> qr((projector_rescale * P_LHP).eval());
    // Eigen's ColPivHouseholderQR::rank() compares each pivot against
    // threshold() times the largest pivot, so setThreshold takes a relative,
    // dimensionless multiplier. The backward-stable rank tolerance for a QR is
    // the matrix size times unit roundoff (2n times epsilon); multiplying by an
    // operand norm would apply the scale twice and inflate the cutoff on
    // exactly the ill-conditioned stable subspaces this path must resolve.
    qr.setThreshold(Scalar{n2} * eps);
    if (qr.rank() < n)
        return ctrlpp::unexpected(care_error::non_lhp_stabilizable);

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

    if (!care_solution_satisfies_postconditions<Scalar, NX>(H_in, out.P))
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
