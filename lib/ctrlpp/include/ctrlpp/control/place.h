#ifndef HPP_GUARD_CTRLPP_CONTROL_PLACE_H
#define HPP_GUARD_CTRLPP_CONTROL_PLACE_H

/// @brief Pole placement via Ackermann's formula for single-input systems.
///
/// @cite kautsky1985 -- Kautsky, Nichols & Van Dooren, "Robust Pole Assignment in Linear State Feedback", 1985

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <optional>

namespace ctrlpp
{

namespace detail
{

/// Validate that complex poles come in conjugate pairs (required for real coefficients).
///
/// A real-coefficient polynomial has roots that are either real or appear in
/// complex-conjugate pairs; any other configuration would force complex
/// coefficients in the closed-loop characteristic polynomial.
///
/// @cite franklin2015 -- Franklin, Powell &amp; Emami-Naeini, "Feedback Control of Dynamic Systems", 2015, Ch. 7
///
/// The conjugate-pair and real-pole tests are relative: a pole's imaginary part
/// and the residual between a candidate pair are judged against the pole
/// magnitude scaled by the instantiated Scalar's machine epsilon and an exposed
/// dimensionless coefficient. This keeps the test correct at both float and
/// double precision and for large-magnitude poles, where a fixed absolute
/// tolerance would either reject valid float pairs or accept spurious double
/// ones. The coefficient defaults to twice the pole count, reflecting the
/// rounding that accumulates across the N pole values being compared, and is
/// overridable by the caller.
template <typename Scalar, std::size_t N>
bool validate_conjugate_pairs(const std::array<std::complex<Scalar>, N>& poles,
                              Scalar conj_tol_scale = Scalar{2} * Scalar{N})
{
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    // Magnitude floor so a zero-magnitude pole still receives an epsilon-level
    // tolerance instead of a tolerance that collapses to zero; the floor is
    // itself the Scalar's epsilon, so no absolute magic constant is introduced.
    const Scalar mag_floor = eps;

    // Track which poles have been matched
    std::array<bool, N> matched{};

    for(std::size_t i = 0; i < N; ++i)
    {
        if(matched[i])
            continue;

        const Scalar tol_i = conj_tol_scale * eps * std::max(std::abs(poles[i]), mag_floor);
        if(std::abs(poles[i].imag()) <= tol_i)
        {
            // Real pole, no conjugate needed
            matched[i] = true;
            continue;
        }

        // Complex pole: find its conjugate
        bool found = false;
        for(std::size_t j = i + 1; j < N; ++j)
        {
            if(matched[j])
                continue;
            const Scalar scale = std::max(std::abs(poles[i]), std::abs(poles[j]));
            const Scalar tol = conj_tol_scale * eps * std::max(scale, mag_floor);
            if(std::abs(poles[i].real() - poles[j].real()) <= tol && std::abs(poles[i].imag() + poles[j].imag()) <= tol)
            {
                matched[i] = true;
                matched[j] = true;
                found = true;
                break;
            }
        }
        if(!found)
            return false;
    }
    return true;
}

/// Compute characteristic polynomial coefficients from desired poles.
/// Returns coefficients [a_0, a_1, ..., a_{N-1}] of:
/// p(s) = s^N + a_{N-1} s^{N-1} + ... + a_1 s + a_0.
///
/// Implements the standard expansion p(s) = prod_i (s - p_i) by repeated
/// polynomial multiplication, which is the constructive form of Vieta's
/// formulas relating elementary symmetric functions of the roots to the
/// signed coefficients.
///
/// @cite strang2016 -- Strang, "Introduction to Linear Algebra", 2016, Ch. 6 (eigenvalues and characteristic polynomial)
template <typename Scalar, std::size_t N>
std::array<Scalar, N> char_poly_coeffs(const std::array<std::complex<Scalar>, N>& poles)
{
    // Start with p(s) = 1, multiply by (s - p_i) one at a time.
    // coeffs[k] stores coefficient of s^k in the accumulated polynomial.
    // We use complex arithmetic internally and take real parts at the end.
    std::array<std::complex<Scalar>, N + 1> c{};
    c[0] = std::complex<Scalar>{1, 0};

    for(std::size_t i = 0; i < N; ++i)
    {
        // Multiply current poly by (s - poles[i])
        // New poly degree is i+1
        // Process from high to low to avoid overwriting
        for(std::size_t k = i + 1; k > 0; --k)
            c[k] = c[k - 1] - poles[i] * c[k];
        c[0] = -poles[i] * c[0];
    }

    // c[N] = 1 (leading coeff), c[k] = coeff of s^k for k=0..N-1
    // Return [a_0, ..., a_{N-1}]
    std::array<Scalar, N> result;
    for(std::size_t k = 0; k < N; ++k)
        result[k] = c[k].real();
    return result;
}

}

/// Pole placement using Ackermann's formula for single-input systems (NU == 1).
/// Computes K such that eigenvalues of (A - B*K) equal the desired poles.
/// Returns std::nullopt if the system is uncontrollable or NU &gt; 1.
///
/// The single-input formula K = e_n^T * C_ctrl^{-1} * alpha(A) follows from
/// transforming (A, B) into controller canonical form, assigning the desired
/// characteristic polynomial alpha(s), and back-transforming the gain.
///
/// @cite kautsky1985 -- Kautsky, Nichols &amp; Van Dooren, "Robust Pole Assignment in Linear State Feedback", 1985
/// @cite franklin2015 -- Franklin, Powell &amp; Emami-Naeini, "Feedback Control of Dynamic Systems", 2015, Ch. 7 (Ackermann's formula)
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>> place(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const std::array<std::complex<Scalar>, NX>& desired_poles, Scalar conj_tol_scale = Scalar{2} * Scalar{NX})
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    constexpr int n = static_cast<int>(NX);

    if constexpr(NU != 1)
    {
        // Multi-input placement not yet supported
        return std::nullopt;
    }
    else
    {
        // Validate conjugate pairs
        if(!detail::validate_conjugate_pairs(desired_poles, conj_tol_scale))
            return std::nullopt;

        // Build controllability matrix: C_ctrl = [B, AB, A^2 B, ..., A^{n-1} B]
        Eigen::Matrix<Scalar, n, n> C_ctrl;
        Eigen::Matrix<Scalar, n, 1> AkB = B;
        C_ctrl.col(0) = AkB;
        for(int k = 1; k < n; ++k)
        {
            AkB = (A * AkB).eval();
            C_ctrl.col(k) = AkB;
        }

        // Check controllability
        auto qr = C_ctrl.colPivHouseholderQr();
        if(qr.rank() < n)
            return std::nullopt;

        // Compute characteristic polynomial coefficients
        auto alpha = detail::char_poly_coeffs(desired_poles);

        // Evaluate alpha(A) using Horner's method:
        // alpha(A) = A^n + alpha_{n-1} A^{n-1} + ... + alpha_1 A + alpha_0 I
        Eigen::Matrix<Scalar, n, n> alphaA = Eigen::Matrix<Scalar, n, n>::Identity();
        // Horner: start from highest power
        for(int k = n - 1; k >= 0; --k)
            alphaA = (alphaA * A + alpha[static_cast<std::size_t>(k)] * Eigen::Matrix<Scalar, n, n>::Identity()).eval();

        // K = e_n^T * C_ctrl^{-1} * alpha(A)
        // e_n^T = last row of identity = [0, 0, ..., 0, 1]
        Eigen::Matrix<Scalar, n, n> C_ctrl_inv = qr.inverse();
        Eigen::Matrix<Scalar, 1, n> en_T = Eigen::Matrix<Scalar, 1, n>::Zero();
        en_T(0, n - 1) = Scalar{1};

        Eigen::Matrix<Scalar, 1, n> K = en_T * C_ctrl_inv * alphaA;
        return K;
    }
}

/// Convenience: compute observer gain L via duality.
/// L = place(A^T, C^T, desired_poles)^T.
///
/// @cite franklin2015 -- Franklin, Powell &amp; Emami-Naeini, "Feedback Control of Dynamic Systems", 2015, Ch. 7 (observer/regulator duality)
template <typename Scalar, std::size_t NX, std::size_t NY>
std::optional<Eigen::Matrix<Scalar, int(NX), int(NY)>>
place_observer(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NY), int(NX)>& C, const std::array<std::complex<Scalar>, NX>& desired_poles, Scalar conj_tol_scale = Scalar{2} * Scalar{NX})
{
    if constexpr(NY != 1)
    {
        // Multi-output observer placement not yet supported (duality requires single-output)
        return std::nullopt;
    }
    else
    {
        Eigen::Matrix<Scalar, int(NX), int(NX)> At = A.transpose().eval();
        Eigen::Matrix<Scalar, int(NX), int(NY)> Ct = C.transpose().eval();

        auto K_opt = place<Scalar, NX, NY>(At, Ct, desired_poles, conj_tol_scale);
        if(!K_opt)
            return std::nullopt;

        return K_opt->transpose().eval();
    }
}

}

#endif
