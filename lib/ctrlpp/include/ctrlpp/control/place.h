#ifndef HPP_GUARD_CTRLPP_CONTROL_PLACE_H
#define HPP_GUARD_CTRLPP_CONTROL_PLACE_H

// Ackermann's formula K = e_n^T C_ctrl^{-1} alpha(A): transform (A, B) into controller canonical form, assign the
// desired characteristic polynomial alpha(s), and transform the gain back.
// @cite kautsky1985 -- Kautsky, Nichols & Van Dooren, "Robust Pole Assignment in Linear State Feedback", 1985
// @cite franklin2015 -- Franklin, Powell & Emami-Naeini, "Feedback Control of Dynamic Systems", 2015, Ch. 7 (Ackermann's formula)

#include "ctrlpp/control/place_types.h"

#include "ctrlpp/control/detail/place_poles.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU>
bool place_inputs_finite(
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const std::array<std::complex<Scalar>, NX>& poles, Scalar conj_tol_scale)
{
    const auto finite_pole = [](const std::complex<Scalar>& pole) { return std::isfinite(pole.real()) && std::isfinite(pole.imag()); };
    return A.allFinite() && B.allFinite() && std::isfinite(conj_tol_scale) && std::all_of(poles.begin(), poles.end(), finite_pole);
}

template <typename Scalar, std::size_t NX>
Eigen::Matrix<Scalar, int(NX), int(NX)> controllability_matrix(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), 1>& B)
{
    Eigen::Matrix<Scalar, int(NX), int(NX)> C_ctrl;
    Eigen::Matrix<Scalar, int(NX), 1> AkB = B;
    C_ctrl.col(0) = AkB;
    for(Eigen::Index k = 1; k < Eigen::Index(NX); ++k)
    {
        AkB = (A * AkB).eval();
        C_ctrl.col(k) = AkB;
    }
    return C_ctrl;
}

// alpha(A) = A^n + alpha_{n-1} A^{n-1} + ... + alpha_0 I by Horner's method.
template <typename Scalar, std::size_t NX>
Eigen::Matrix<Scalar, int(NX), int(NX)> characteristic_polynomial_of(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const std::array<std::complex<Scalar>, NX>& poles)
{
    using MatN = Eigen::Matrix<Scalar, int(NX), int(NX)>;

    const std::array<Scalar, NX> alpha = char_poly_coeffs(poles);
    MatN alphaA = MatN::Identity();
    for(std::size_t k = NX; k > 0; --k)
        alphaA = (alphaA * A + alpha[k - 1] * MatN::Identity()).eval();
    return alphaA;
}

template <typename Scalar, std::size_t NX>
Eigen::Matrix<Scalar, 1, int(NX)> ackermann_gain(const Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, int(NX), int(NX)>>& qr, const Eigen::Matrix<Scalar, int(NX), int(NX)>& alphaA)
{
    constexpr int n = static_cast<int>(NX);

    Eigen::Matrix<Scalar, n, n> C_ctrl_inv = qr.inverse();
    Eigen::Matrix<Scalar, 1, n> en_T = Eigen::Matrix<Scalar, 1, n>::Zero();
    en_T(0, n - 1) = Scalar{1};

    Eigen::Matrix<Scalar, 1, n> K = en_T * C_ctrl_inv * alphaA;
    return K;
}

}

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, place_error> place(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const std::array<std::complex<Scalar>, NX>& desired_poles, Scalar conj_tol_scale = Scalar{2} * Scalar{NX})
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    if constexpr(NU != 1)
        return ctrlpp::unexpected(place_error::multi_input_not_supported);
    else
    {
        if(!detail::place_inputs_finite<Scalar, NX, NU>(A, B, desired_poles, conj_tol_scale))
            return ctrlpp::unexpected(place_error::non_finite_input);
        if(!detail::validate_conjugate_pairs(desired_poles, conj_tol_scale))
            return ctrlpp::unexpected(place_error::poles_not_conjugate_symmetric);

        const Eigen::Matrix<Scalar, int(NX), int(NX)> C_ctrl = detail::controllability_matrix<Scalar, NX>(A, B);
        auto qr = C_ctrl.colPivHouseholderQr();
        if(qr.rank() < Eigen::Index(NX))
            return ctrlpp::unexpected(place_error::uncontrollable_pair);

        return detail::ackermann_gain<Scalar, NX>(qr, detail::characteristic_polynomial_of<Scalar, NX>(A, desired_poles));
    }
}

// @cite franklin2015 -- Franklin, Powell & Emami-Naeini, "Feedback Control of Dynamic Systems", 2015, Ch. 7 (observer/regulator duality)
template <typename Scalar, std::size_t NX, std::size_t NY>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NX), int(NY)>, place_error> place_observer(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NY), int(NX)>& C,
        const std::array<std::complex<Scalar>, NX>& desired_poles, Scalar conj_tol_scale = Scalar{2} * Scalar{NX})
{
    if constexpr(NY != 1)
        return ctrlpp::unexpected(place_error::multi_output_not_supported);
    else
    {
        Eigen::Matrix<Scalar, int(NX), int(NX)> At = A.transpose().eval();
        Eigen::Matrix<Scalar, int(NX), int(NY)> Ct = C.transpose().eval();

        auto K_opt = place<Scalar, NX, NY>(At, Ct, desired_poles, conj_tol_scale);
        if(!K_opt)
            return ctrlpp::unexpected(K_opt.error());

        return Eigen::Matrix<Scalar, int(NX), int(NY)>{K_opt->transpose().eval()};
    }
}

}

#endif
