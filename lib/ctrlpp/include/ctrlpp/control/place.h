#ifndef HPP_GUARD_CTRLPP_CONTROL_PLACE_H
#define HPP_GUARD_CTRLPP_CONTROL_PLACE_H

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>

namespace ctrlpp {

namespace detail {

// Validate that complex poles come in conjugate pairs (required for real coefficients).
template <typename Scalar, std::size_t N>
bool validate_conjugate_pairs(const std::array<std::complex<Scalar>, N> &poles)
{
    // Track which poles have been matched
    std::array<bool, N> matched{};

    for(std::size_t i = 0; i < N; ++i)
    {
        if(matched[i])
            continue;

        if(std::abs(poles[i].imag()) < Scalar{1e-12})
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
            if(std::abs(poles[i].real() - poles[j].real()) < Scalar{1e-12} && std::abs(poles[i].imag() + poles[j].imag()) < Scalar{1e-12})
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

// Compute characteristic polynomial coefficients from desired poles.
// Returns coefficients [a_0, a_1, ..., a_{N-1}] of:
// p(s) = s^N + a_{N-1} s^{N-1} + ... + a_1 s + a_0
template <typename Scalar, std::size_t N>
std::array<Scalar, N> char_poly_coeffs(const std::array<std::complex<Scalar>, N> &poles)
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

// Pole placement using Ackermann's formula for single-input systems (NU == 1).
// Computes K such that eigenvalues of (A - B*K) equal the desired poles.
// Returns std::nullopt if the system is uncontrollable or NU > 1.
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>> place(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B, const std::array<std::complex<Scalar>, NX> &desired_poles)
{
    constexpr int n = static_cast<int>(NX);

    if constexpr(NU != 1)
    {
        // Multi-input placement not yet supported
        return std::nullopt;
    }
    else
    {
        // Validate conjugate pairs
        if(!detail::validate_conjugate_pairs(desired_poles))
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
            alphaA = (alphaA * A + alpha[static_cast<std::size_t>(k)] *
                Eigen::Matrix<Scalar, n, n>::Identity()).eval();

        // K = e_n^T * C_ctrl^{-1} * alpha(A)
        // e_n^T = last row of identity = [0, 0, ..., 0, 1]
        Eigen::Matrix<Scalar, n, n> C_ctrl_inv = qr.inverse();
        Eigen::Matrix<Scalar, 1, n> en_T = Eigen::Matrix<Scalar, 1, n>::Zero();
        en_T(0, n - 1) = Scalar{1};

        Eigen::Matrix<Scalar, 1, n> K = en_T * C_ctrl_inv * alphaA;
        return K;
    }
}

// Convenience: compute observer gain L via duality.
// L = place(A^T, C^T, desired_poles)^T
template <typename Scalar, std::size_t NX, std::size_t NY>
std::optional<Eigen::Matrix<Scalar, int(NX), int(NY)>> place_observer(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NY), int(NX)> &C, const std::array<std::complex<Scalar>, NX> &desired_poles)
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

        auto K_opt = place<Scalar, NX, NY>(At, Ct, desired_poles);
        if(!K_opt)
            return std::nullopt;

        return K_opt->transpose().eval();
    }
}

}

#endif
