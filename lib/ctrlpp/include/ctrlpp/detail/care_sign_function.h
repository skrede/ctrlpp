#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_SIGN_FUNCTION_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_SIGN_FUNCTION_H

// Continuous Riccati solve by the scaled Newton iteration H_{k+1} = (mu_k H_k + mu_k^-1 H_k^-1) / 2 on the Hamiltonian,
// which converges to sign(H). The determinantal scaling mu_k = |det H_k|^(-1/(2n)) is formed overflow-safe from the LU
// diagonal as exp(-(1/(2n)) sum_i log|u_ii|). The stable invariant subspace is the range of (I - sign(H)) / 2, whose
// rank-revealing QR gives the basis the shared extraction reads. A small change between iterates is only a candidate for
// convergence: the shared postconditions on the returned solution decide acceptance, and the projector's idempotence is
// not consulted, because it holds exactly even where the extracted solution comes back as zero.
// @cite roberts1980 : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
// @cite byers1987   : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
// @cite higham2008  : Higham, "Functions of Matrices: Theory and Computation", 2008, Ch. 5 Eq. 5.16 / 5.34 / 5.35 / 5.41, Theorem 5.1(e)

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/care_postconditions.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/control/care_types.h"

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

// The change test decides only when to stop. The iteration converges quadratically, so a step of relative change d leaves
// an error of order d^2; setting d^2 to the counted floor of the 2n x 2n difference's Frobenius magnitude, n2^2 eps,
// gives d = n2 sqrt(eps).
template <typename Scalar, int N2>
Scalar sign_change_tolerance()
{
    constexpr std::int32_t change_candidate_rounding_ops = N2 * N2;
    return std::sqrt(static_cast<Scalar>(change_candidate_rounding_ops) * std::numeric_limits<Scalar>::epsilon());
}

template <typename Scalar, int N2>
ctrlpp::expected<void, care_error> scaled_newton_step(Eigen::Matrix<Scalar, N2, N2>& H)
{
    using Mat2N = Eigen::Matrix<Scalar, N2, N2>;
    Eigen::PartialPivLU<Mat2N> lu(H);
    Scalar log_det_abs = Scalar{0};
    for(Eigen::Index i = 0; i < N2; ++i)
        log_det_abs += std::log(std::abs(lu.matrixLU()(i, i)));
    const Scalar mu = std::exp(-log_det_abs / Scalar{N2});

    if(!std::isfinite(mu))
        return ctrlpp::unexpected(care_error::sign_function_stagnated);

    const Mat2N H_inv = lu.solve(Mat2N::Identity()).eval();
    H = ((Scalar{1} / Scalar{2}) * (mu * H + (Scalar{1} / mu) * H_inv)).eval();

    if(!H.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);
    return {};
}

template <typename Scalar, int N2>
ctrlpp::expected<bool, care_error> sign_iteration_converged(
        const Eigen::Matrix<Scalar, N2, N2>& H, const Eigen::Matrix<Scalar, N2, N2>& H_prev, Scalar tolerance, bool guard_armed, resolved_magnitude<Scalar>& last_change)
{
    // On the first step the previous iterate is the caller's Hamiltonian, whose sum of squares can overflow while every
    // entry is finite, so both operands are resolved and an unresolved one leaves through the typed channel.
    const resolved_magnitude<Scalar> change = resolve_magnitude(H - H_prev);
    const resolved_magnitude<Scalar> scale_H = resolve_magnitude(H);
    if(!change.resolved || !scale_H.resolved)
        return ctrlpp::unexpected(care_error::sign_function_stagnated);

    const resolved_magnitude<Scalar> change_candidate_floor = scaled_magnitude(scale_H, tolerance);
    if(magnitude_within(change, change_candidate_floor))
        return true;

    // A change that grows means the iteration is not contracting toward a sign matrix, so the solve declines now rather
    // than spending the rest of the budget. The determinantal scaling makes large early corrections, which is why the
    // guard waits out the warmup window.
    if(guard_armed && !magnitude_within(change, last_change))
        return ctrlpp::unexpected(care_error::sign_function_stagnated);
    last_change = change;
    return false;
}

template <typename Scalar, int N2>
ctrlpp::expected<void, care_error> sign_function_projector(const Eigen::Matrix<Scalar, N2, N2>& H_in, std::int32_t warmup_iterations, Eigen::Matrix<Scalar, N2, N2>& P_LHP)
{
    using Mat2N = Eigen::Matrix<Scalar, N2, N2>;
    // Higham Table 5.2's worst observed iteration count with determinantal scaling.
    constexpr std::int32_t max_iters = 40;
    const Scalar change_candidate_tolerance = sign_change_tolerance<Scalar, N2>();
    resolved_magnitude<Scalar> last_change{Scalar{0}, false};
    Mat2N H = H_in;

    for(std::int32_t k = 0; k < max_iters; ++k)
    {
        const Mat2N H_prev = H;
        const ctrlpp::expected<void, care_error> step = scaled_newton_step<Scalar, N2>(H);
        if(!step)
            return step;
        const ctrlpp::expected<bool, care_error> converged = sign_iteration_converged<Scalar, N2>(H, H_prev, change_candidate_tolerance, k > warmup_iterations, last_change);
        if(!converged)
            return ctrlpp::unexpected(converged.error());
        if(*converged)
        {
            P_LHP = ((Scalar{1} / Scalar{2}) * (Mat2N::Identity() - H)).eval();
            return {};
        }
    }
    return ctrlpp::unexpected(care_error::sign_function_stagnated);
}

// Each Householder step drops a reflection whose column tail has a squared norm at or below the smallest normal value,
// which zeroed solutions below its square root. The QR is equivariant under a positive scale, so the projector is scaled
// by 2^(e - 2), where frexp writes sqrt(max / n2) / max|P_LHP| as m 2^e: an exact power of two that bounds every tail's
// sum of squares strictly below a quarter of the largest finite value.
template <typename Scalar, int N2>
Scalar projector_rescale_factor(const Eigen::Matrix<Scalar, N2, N2>& P_LHP)
{
    Scalar projector_rescale = Scalar{1};
    const resolved_magnitude<Scalar> projector_largest_entry = resolve_largest_entry(P_LHP);
    if(projector_largest_entry.resolved && projector_largest_entry.value > Scalar{0})
    {
        const Scalar tail_ceiling = std::sqrt(std::numeric_limits<Scalar>::max() / Scalar{N2});
        const Scalar headroom_ratio = tail_ceiling / projector_largest_entry.value;
        if(std::isfinite(headroom_ratio) && headroom_ratio > Scalar{0})
        {
            int headroom_exponent = 0;
            std::frexp(headroom_ratio, &headroom_exponent);
            const Scalar candidate = std::ldexp(Scalar{1}, headroom_exponent - 2);
            if(std::isfinite(candidate) && candidate > Scalar{0})
                projector_rescale = candidate;
        }
    }
    return projector_rescale;
}

template <typename Scalar, std::size_t NX>
ctrlpp::expected<void, care_error> extract_sign_function_solution(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& P_LHP, Eigen::Matrix<Scalar, int(NX), int(NX)>& P)
{
    constexpr int n = int(NX);
    constexpr int n2 = 2 * n;

    const Scalar projector_rescale = projector_rescale_factor<Scalar, n2>(P_LHP);
    Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, n2, n2>> qr((projector_rescale * P_LHP).eval());
    // rank() compares each pivot against threshold() times the largest pivot, so the threshold is the relative multiplier
    // 2n eps; multiplying it by an operand norm would apply the scale twice.
    qr.setThreshold(Scalar{n2} * std::numeric_limits<Scalar>::epsilon());
    if(qr.rank() < n)
        return ctrlpp::unexpected(care_error::non_lhp_stabilizable);

    const Eigen::Matrix<Scalar, n2, n2> U = qr.householderQ();
    const ctrlpp::expected<void, riccati_extract_error> P_err = extract_riccati_solution_into<Scalar, n2>(P, U);
    if(!P_err)
        return ctrlpp::unexpected(care_error_from_extract(P_err.error()));
    return {};
}

template <typename Scalar, std::size_t NX>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care_solve_via_sign_function(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in, std::int32_t warmup_iterations)
{
    constexpr int n2 = 2 * int(NX);

    if(!H_in.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    Eigen::Matrix<Scalar, n2, n2> P_LHP;
    const ctrlpp::expected<void, care_error> projector = sign_function_projector<Scalar, n2>(H_in, warmup_iterations, P_LHP);
    if(!projector)
        return ctrlpp::unexpected(projector.error());

    care_result<Scalar, NX> out;
    const ctrlpp::expected<void, care_error> extracted = extract_sign_function_solution<Scalar, NX>(P_LHP, out.P);
    if(!extracted)
        return ctrlpp::unexpected(extracted.error());
    if(!care_solution_satisfies_postconditions<Scalar, NX>(H_in, out.P))
        return ctrlpp::unexpected(care_error::sign_function_stagnated);

    // No swap phase, so no pivot ratio: quiet_NaN marks it unavailable, as the care_result contract states.
    out.subspace_separation = std::numeric_limits<Scalar>::quiet_NaN();
    out.reorder_complete = true;
    return out;
}

}

#endif
