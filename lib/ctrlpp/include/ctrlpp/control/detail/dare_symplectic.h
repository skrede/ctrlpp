#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_DARE_SYMPLECTIC_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_DARE_SYMPLECTIC_H

#include "ctrlpp/control/dare_types.h"

#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/power_of_two_scale.h"

#include "ctrlpp/expected.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <limits>
#include <complex>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
struct dare_symplectic_operands
{
    Eigen::Matrix<Scalar, int(NX), int(NX)> A_inverse_transpose;
    Eigen::Matrix<Scalar, int(NX), int(NX)> G;
};

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NX), int(NX)>, dare_error> factor_dare_g(const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    // A rank-deficient R would make solve() return a least-squares result and silently pose a different control
    // problem, so the reciprocal-pivot test comes first.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_R.isInvertible())
        return ctrlpp::unexpected(dare_error::singular_r);

    return Eigen::Matrix<Scalar, int(NX), int(NX)>{(B * qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval()};
}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<dare_symplectic_operands<Scalar, NX>, dare_error> factor_dare_symplectic_operands(
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    constexpr int n = static_cast<int>(NX);
    using MatNxN = Eigen::Matrix<Scalar, n, n>;

    auto qr_At = A.transpose().colPivHouseholderQr();
    qr_At.setThreshold(Scalar{static_cast<int>(NX)} * std::numeric_limits<Scalar>::epsilon());
    if(!qr_At.isInvertible())
        return ctrlpp::unexpected(dare_error::singular_a);

    auto G_result = factor_dare_g<Scalar, NX, NU>(B, R);
    if(!G_result)
        return ctrlpp::unexpected(G_result.error());

    dare_symplectic_operands<Scalar, NX> operands;
    operands.A_inverse_transpose = qr_At.solve(MatNxN::Identity()).eval();
    operands.G = *G_result;
    return operands;
}

// A common positive divisor of Q and R preserves the gain, and a power of two makes the division of the weights and the
// multiplication of the solution back exact.
template <typename Scalar, std::size_t NX, std::size_t NU>
Scalar dare_weight_scale(const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    return power_of_two_ceiling(std::max(Q.cwiseAbs().maxCoeff(), R.cwiseAbs().maxCoeff()));
}

// The symplectic matrix of Laub 1979 Eq. 7, with G = B R^{-1} B'.
template <typename Scalar, std::size_t NX>
ctrlpp::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, dare_error> build_dare_symplectic(
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const dare_symplectic_operands<Scalar, NX>& operands)
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;

    Mat2Nx2N Z;
    Z.template block<n, n>(0, 0) = A + operands.G * operands.A_inverse_transpose * Q;
    Z.template block<n, n>(0, n) = -operands.G * operands.A_inverse_transpose;
    Z.template block<n, n>(n, 0) = -operands.A_inverse_transpose * Q;
    Z.template block<n, n>(n, n) = operands.A_inverse_transpose;

    if(!Z.allFinite())
        return ctrlpp::unexpected(dare_error::arithmetic_limit);

    return Z;
}

template <typename Scalar, std::size_t NX, conditioning_policy Cond>
ctrlpp::expected<reorder_result<Scalar>, dare_error> reorder_dare_schur(
        const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& Z, Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& T, Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& U)
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;

    Eigen::RealSchur<Eigen::Matrix<Scalar, n2, n2>> schur(Z);
    if(schur.info() != Eigen::Success)
        return ctrlpp::unexpected(dare_error::schur_failed);

    T = schur.matrixT();
    U = schur.matrixU();
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);

    const Scalar scale = T.cwiseAbs().maxCoeff();
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    // The backward error of the real Schur factor's eigenvalues, 2n eps ||T||_max, is the unit-disk margin.
    const Scalar unit_margin = Scalar{n2} * eps * scale;
    auto predicate = [unit_margin](std::complex<Scalar> lam) -> bool { return std::abs(lam) < Scalar{1} - unit_margin; };

    auto rr = reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
    if(rr.placed < n)
        return ctrlpp::unexpected(dare_error::non_stabilizable);
    if(!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(dare_error::non_finite_input);
    return rr;
}

inline dare_error dare_error_from_extract(riccati_extract_error error)
{
    switch(error)
    {
        case riccati_extract_error::singular_u11:
            return dare_error::singular_u11;
        case riccati_extract_error::non_finite:
            return dare_error::non_finite_input;
        case riccati_extract_error::non_psd:
            return dare_error::non_psd_solution;
    }
    return dare_error::non_finite_input;
}

// The result's gain is left zero: this entry point sees only Z, and dare fills the gain from the verification that
// decided the acceptance.
template <typename Scalar, std::size_t NX, std::size_t NU, conditioning_policy Cond = pivot_ratio_conditioning>
ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error> dare_solve_from_symplectic(const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& Z, Cond /*tag*/ = {})
{
    constexpr int n2 = 2 * static_cast<int>(NX);
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

    Mat2N T;
    Mat2N U;
    const ctrlpp::expected<reorder_result<Scalar>, dare_error> rr = reorder_dare_schur<Scalar, NX, Cond>(Z, T, U);
    if(!rr)
        return ctrlpp::unexpected(rr.error());

    dare_result<Scalar, NX, NU> out;
    auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
    if(!P_err)
        return ctrlpp::unexpected(dare_error_from_extract(P_err.error()));

    out.K.setZero();
    out.subspace_separation = rr->subspace_separation;
    out.reorder_complete = rr->complete;
    return out;
}

}

#endif
