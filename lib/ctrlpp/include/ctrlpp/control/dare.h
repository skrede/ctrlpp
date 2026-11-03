#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_H

// Schur method on the symplectic matrix, with the eigenvalues inside the unit disk reordered to the leading block by the
// Bai-Demmel swap kernel and P = U21 U11^-1 extracted from the resulting invariant subspace basis.
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include "ctrlpp/control/dare_types.h"

#include "ctrlpp/control/detail/dare_symplectic.h"
#include "ctrlpp/control/detail/dare_verification.h"

#include "ctrlpp/detail/schur_reorder.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU>
bool equilibrate_dare_weights(Scalar weight_scale, Eigen::Matrix<Scalar, int(NX), int(NX)>& Q_scaled, Eigen::Matrix<Scalar, int(NU), int(NU)>& R_scaled)
{
    if(!(weight_scale > Scalar{0}) || !std::isfinite(weight_scale))
        return false;
    if(weight_scale == Scalar{1})
        return true;

    Q_scaled /= weight_scale;
    R_scaled /= weight_scale;
    return Q_scaled.allFinite() && R_scaled.allFinite();
}

// Equilibration exists to make the check's quantities formable, so at this scale only a verified verdict is accepted
// and gain_unavailable declines. The gain is homogeneous of degree zero in (P, Q, R), so the equilibrated gain is the
// caller's gain unscaled.
template <typename Scalar, std::size_t NX, std::size_t NU, conditioning_policy Cond>
ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error> solve_equilibrated_dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q_scaled, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R_scaled)
{
    auto operands_result = factor_dare_symplectic_operands<Scalar, NX, NU>(A, B, R_scaled);
    if(!operands_result)
        return ctrlpp::unexpected(operands_result.error());

    auto Z_result = build_dare_symplectic<Scalar, NX>(A, Q_scaled, *operands_result);
    if(!Z_result)
        return ctrlpp::unexpected(Z_result.error());

    auto result = dare_solve_from_symplectic<Scalar, NX, NU, Cond>(*Z_result);
    if(!result)
        return ctrlpp::unexpected(result.error() == dare_error::non_finite_input ? dare_error::arithmetic_limit : result.error());

    if(!Q_scaled.allFinite() || !R_scaled.allFinite() || verify_dare_solution<Scalar, NX, NU>(A, B, Q_scaled, R_scaled, result->P, result->K) != dare_verification::verified)
        return ctrlpp::unexpected(dare_error::arithmetic_limit);
    return result;
}

// The reduction inverts R before the symplectic build sees it; without this test a singular R would reach the inner
// solve as a non-finite Q' and A' and be reported as a non-finite input.
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error> dare_cross_weight_correction(
        const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NU)>& N)
{
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_R.isInvertible())
        return ctrlpp::unexpected(dare_error::singular_r);

    return Eigen::Matrix<Scalar, int(NU), int(NX)>{qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose())).eval()};
}

}

// Equilibrates before factorizing: on a subnormal weighting R^{-1} overflows and B R^{-1} B' is infinite at the caller's
// scale while finite at the equilibrated one, so only the equilibrated R is ever factorized.
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error> dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, Cond /*tag*/ = {})
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    const Scalar weight_scale = detail::dare_weight_scale<Scalar, NX, NU>(Q, R);
    auto Q_scaled = Q;
    auto R_scaled = R;
    if(!detail::equilibrate_dare_weights<Scalar, NX, NU>(weight_scale, Q_scaled, R_scaled))
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    auto result = detail::solve_equilibrated_dare<Scalar, NX, NU, Cond>(A, B, Q_scaled, R_scaled);
    if(!result || weight_scale == Scalar{1})
        return result;

    if(!detail::carry_dare_solution<Scalar, NX, NU>(A, B, Q, R, weight_scale, *result))
        return ctrlpp::unexpected(dare_error::arithmetic_limit);
    return result;
}

// Reduces to standard form with Q' = Q - N R^{-1} N' and A' = A - B R^{-1} N'. P is shared by both posings, and the
// posed problem's gain is K = K' + R^{-1} N'.
template <typename Scalar, std::size_t NX, std::size_t NU, detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error> dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NU)>& N, Cond tag = {})
{
    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite() || !N.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    const auto Rinv_Nt = detail::dare_cross_weight_correction<Scalar, NX, NU>(R, N);
    if(!Rinv_Nt)
        return ctrlpp::unexpected(Rinv_Nt.error());

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * *Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * *Rinv_Nt).eval();
    if(!Qp.allFinite() || !Ap.allFinite())
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    auto reduced = dare<Scalar, NX, NU, Cond>(Ap, B, Qp, R, tag);
    if(!reduced)
        return reduced;

    reduced->K += *Rinv_Nt;
    if(!reduced->K.allFinite())
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    return reduced;
}

}

#endif
