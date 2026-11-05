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

template <typename Scalar, std::size_t NX>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care_solve_via_sign_function(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in, std::int32_t warmup_iterations)
{
    constexpr int n = int(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    if(!H_in.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    Mat2N H = H_in;
    Mat2N H_prev;
    Mat2N P_LHP;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    // Higham Table 5.2's worst observed iteration count with determinantal scaling.
    constexpr std::int32_t max_iters = 40;
    // The change test decides only when to stop. The iteration converges quadratically, so a step of relative change d
    // leaves an error of order d^2; setting d^2 to the counted floor of the 2n x 2n difference's Frobenius magnitude,
    // n2^2 eps, gives d = n2 sqrt(eps).
    constexpr std::int32_t change_candidate_rounding_ops = n2 * n2;
    const Scalar change_candidate_tolerance = std::sqrt(static_cast<Scalar>(change_candidate_rounding_ops) * eps);
    resolved_magnitude<Scalar> last_change{Scalar{0}, false};

    for(std::int32_t k = 0; k < max_iters; ++k)
    {
        H_prev = H;

        Eigen::PartialPivLU<Mat2N> lu(H);
        Scalar log_det_abs = Scalar{0};
        for(Eigen::Index i = 0; i < n2; ++i)
            log_det_abs += std::log(std::abs(lu.matrixLU()(i, i)));
        const Scalar mu = std::exp(-log_det_abs / Scalar{n2});

        if(!std::isfinite(mu))
            return ctrlpp::unexpected(care_error::sign_function_stagnated);

        const Mat2N H_inv = lu.solve(Mat2N::Identity()).eval();
        H = ((Scalar{1} / Scalar{2}) * (mu * H + (Scalar{1} / mu) * H_inv)).eval();

        if(!H.allFinite())
            return ctrlpp::unexpected(care_error::non_finite_input);

        // On the first step the previous iterate is the caller's Hamiltonian, whose sum of squares can overflow while
        // every entry is finite, so both operands are resolved and an unresolved one leaves through the typed channel.
        const resolved_magnitude<Scalar> change = resolve_magnitude(H - H_prev);
        const resolved_magnitude<Scalar> scale_H = resolve_magnitude(H);
        if(!change.resolved || !scale_H.resolved)
            return ctrlpp::unexpected(care_error::sign_function_stagnated);

        const resolved_magnitude<Scalar> change_candidate_floor = scaled_magnitude(scale_H, change_candidate_tolerance);
        if(magnitude_within(change, change_candidate_floor))
        {
            P_LHP = ((Scalar{1} / Scalar{2}) * (Mat2N::Identity() - H)).eval();
            break;
        }

        // A change that grows means the iteration is not contracting toward a sign matrix, so the solve declines now
        // rather than spending the rest of the budget. The determinantal scaling makes large early corrections, which is
        // why the guard waits out the warmup window.
        if(k > warmup_iterations && !magnitude_within(change, last_change))
            return ctrlpp::unexpected(care_error::sign_function_stagnated);
        last_change = change;

        if(k + 1 == max_iters)
            return ctrlpp::unexpected(care_error::sign_function_stagnated);
    }

    // Each Householder step drops a reflection whose column tail has a squared norm at or below the smallest normal
    // value, which zeroed solutions below its square root. The QR is equivariant under a positive scale, so the
    // projector is scaled by 2^(e - 2), where frexp writes sqrt(max / n2) / max|P_LHP| as m 2^e: an exact power of two
    // that bounds every tail's sum of squares strictly below a quarter of the largest finite value.
    Scalar projector_rescale = Scalar{1};
    const resolved_magnitude<Scalar> projector_largest_entry = resolve_largest_entry(P_LHP);
    if(projector_largest_entry.resolved && projector_largest_entry.value > Scalar{0})
    {
        const Scalar tail_ceiling = std::sqrt(std::numeric_limits<Scalar>::max() / Scalar{n2});
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

    Eigen::ColPivHouseholderQR<Mat2N> qr((projector_rescale * P_LHP).eval());
    // rank() compares each pivot against threshold() times the largest pivot, so the threshold is the relative multiplier
    // 2n eps; multiplying it by an operand norm would apply the scale twice.
    qr.setThreshold(Scalar{n2} * eps);
    if(qr.rank() < n)
        return ctrlpp::unexpected(care_error::non_lhp_stabilizable);

    const Mat2N Q_full = qr.householderQ();
    const Mat2N U = Q_full;

    care_result<Scalar, NX> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if(!P_err)
    {
        switch(P_err.error())
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

    if(!care_solution_satisfies_postconditions<Scalar, NX>(H_in, out.P))
        return ctrlpp::unexpected(care_error::sign_function_stagnated);

    // No swap phase, so no pivot ratio: quiet_NaN marks it unavailable, as the care_result contract states.
    out.subspace_separation = std::numeric_limits<Scalar>::quiet_NaN();
    out.reorder_complete = true;
    return out;
}

}

#endif
