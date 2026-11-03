#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_H

// The stabilizing P of A^T P + P A - P B R^{-1} B^T P + Q = 0 from the Hamiltonian H = [[A, -B R^{-1} B^T], [-Q, -A^T]]:
// by default a scaled Newton iteration on sign(H), otherwise real Schur with Bai-Demmel reordering, optionally balanced;
// every method verifies the extracted P against H before reporting success.
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/control/detail/care_hamiltonian.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <limits>
#include <cstddef>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, care_solve_method Method, conditioning_policy Cond>
ctrlpp::expected<care_result<Scalar, NX>, care_error> solve_equilibrated_care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q_scaled, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R_scaled)
{
    auto H_result = build_care_hamiltonian<Scalar, NX, NU>(A, B, Q_scaled, R_scaled);
    if(!H_result)
        return ctrlpp::unexpected(H_result.error());
    return care_solve_from_hamiltonian<Scalar, NX, Method, Cond>(*H_result);
}

}

// The equation is homogeneous of degree one in (P, Q, R), so equilibrating the weightings changes which poses have an
// answer but not the answer.
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, detail::care_solve_method Method = detail::sign_function_care_method,
        detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, Method /*method_tag*/ = {}, Cond /*cond_tag*/ = {})
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    if(!Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    const Scalar weight_scale = detail::care_weight_scale<Scalar, NX, NU>(Q, R);
    auto Q_scaled = Q;
    auto R_scaled = R;
    if(!detail::equilibrate_care_weights<Scalar, NX, NU>(weight_scale, Q_scaled, R_scaled))
        return ctrlpp::unexpected(care_error::non_finite_input);

    auto result = detail::solve_equilibrated_care<Scalar, NX, NU, Method, Cond>(A, B, Q_scaled, R_scaled);
    if(!result || weight_scale == Scalar{1})
        return result;

    // The postcondition ran on the equilibrated Hamiltonian, which is similar to the caller's; only the
    // multiplication back comes after it, so the product is checked.
    result->P *= weight_scale;
    if(!result->P.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);
    return result;
}

// Reduces to standard form with Q' = Q - N R^{-1} N^T and A' = A - B R^{-1} N^T.
template <typename Scalar, std::size_t NX, std::size_t NU, detail::care_solve_method Method = detail::sign_function_care_method,
        detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NU)>& N, Method method_tag = {},
        Cond cond_tag = {})
{
    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite() || !N.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // The reduction inverts R before the Hamiltonian build sees it, so the build's test is applied here.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_R.isInvertible())
        return ctrlpp::unexpected(care_error::singular_r);

    auto Rinv_Nt = qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose())).eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return care<Scalar, NX, NU, Method, Cond>(Ap, B, Qp, R, method_tag, cond_tag);
}

}

#endif
