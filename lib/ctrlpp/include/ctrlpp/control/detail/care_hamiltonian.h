#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_CARE_HAMILTONIAN_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_CARE_HAMILTONIAN_H

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/care_methods.h"
#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/care_sign_function.h"
#include "ctrlpp/detail/power_of_two_scale.h"
#include "ctrlpp/detail/care_postconditions.h"
#include "ctrlpp/detail/hamiltonian_balance.h"

#include "ctrlpp/expected.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

// With D = diag(I, sI), D^-1 H D is the Hamiltonian with both weightings divided by s, so the rescale is a similarity
// that cannot move the spectrum; a power of two puts both off-diagonal blocks at unit order and round-trips exactly.
template <typename Scalar, std::size_t NX, std::size_t NU>
Scalar care_weight_scale(const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    return power_of_two_ceiling(std::max(Q.cwiseAbs().maxCoeff(), R.cwiseAbs().maxCoeff()));
}

// An equilibration that cannot be computed is an absence of evidence, which declines rather than proceeds.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool equilibrate_care_weights(Scalar weight_scale, Eigen::Matrix<Scalar, int(NX), int(NX)>& Q_scaled, Eigen::Matrix<Scalar, int(NU), int(NU)>& R_scaled)
{
    if(!(weight_scale > Scalar{0}) || !std::isfinite(weight_scale))
        return false;
    if(weight_scale == Scalar{1})
        return true;

    Q_scaled /= weight_scale;
    R_scaled /= weight_scale;
    return Q_scaled.allFinite() && R_scaled.allFinite();
}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, care_error> build_care_hamiltonian(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
        const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    constexpr int n = static_cast<int>(NX);

    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // Without the reciprocal-pivot test a rank-deficient R gets a least-squares solve, and the Hamiltonian then
    // poses a different problem from the one the caller defined.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_R.isInvertible())
        return ctrlpp::unexpected(care_error::singular_r);
    const Eigen::Matrix<Scalar, n, n> S = (B * qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval();

    Eigen::Matrix<Scalar, 2 * n, 2 * n> H;
    H.template block<n, n>(0, 0) = A;
    H.template block<n, n>(0, n) = -S;
    H.template block<n, n>(n, 0) = -Q;
    H.template block<n, n>(n, n) = -A.transpose();
    if(!H.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);
    return H;
}

// Leaves the stable-subspace basis of H in the leading n columns of U.
template <typename Scalar, std::size_t NX, typename Cond>
ctrlpp::expected<reorder_result<Scalar>, care_error> stable_schur_basis(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H, Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& U)
{
    constexpr int n2 = 2 * int(NX);
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    Eigen::RealSchur<Mat2N> schur(H);
    if(schur.info() != Eigen::Success)
        return ctrlpp::unexpected(care_error::schur_failed);

    Mat2N T = schur.matrixT();
    U = schur.matrixU();
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    const reorder_result<Scalar> rr = reorder_left_half_plane<Scalar, n2, Cond>(T, U);
    if(rr.placed < int(NX))
        return ctrlpp::unexpected(care_error::non_lhp_stabilizable);
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);
    return rr;
}

template <typename Scalar, std::size_t NX, typename Cond>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care_solve_via_schur(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H)
{
    constexpr int n2 = 2 * int(NX);

    Eigen::Matrix<Scalar, n2, n2> U;
    const ctrlpp::expected<reorder_result<Scalar>, care_error> rr = stable_schur_basis<Scalar, NX, Cond>(H, U);
    if(!rr)
        return ctrlpp::unexpected(rr.error());

    care_result<Scalar, NX> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if(!P_err)
        return ctrlpp::unexpected(care_error_from_extract(P_err.error()));

    // Verified against H, as the sign-function branch is, not against T or U, which describe an invariant subspace
    // rather than the equation P must satisfy.
    if(!care_solution_satisfies_postconditions<Scalar, NX>(H, out.P))
        return ctrlpp::unexpected(care_error::unverified_solution);

    out.subspace_separation = rr->subspace_separation;
    out.reorder_complete = rr->complete;
    return out;
}

template <typename Scalar, std::size_t NX, care_solve_method Method = sign_function_care_method, conditioning_policy Cond = pivot_ratio_conditioning>
ctrlpp::expected<care_result<Scalar, NX>, care_error> care_solve_from_hamiltonian(
        const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H, Method method_tag = {}, Cond /*cond_tag*/ = {})
{
    if constexpr(std::is_same_v<Method, schur_care_method>)
        return care_solve_via_schur<Scalar, NX, Cond>(H);
    else if constexpr(std::is_same_v<Method, sign_function_care_method>)
        return care_solve_via_sign_function<Scalar, NX>(H, method_tag.warmup_iterations);
    else
        return care_solve_via_balanced_schur<Scalar, NX, Cond>(H);
}

}

#endif
