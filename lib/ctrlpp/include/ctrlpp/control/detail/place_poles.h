#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PLACE_POLES_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PLACE_POLES_H

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

// The tolerance is relative to the pole magnitude, floored at epsilon so a zero pole still gets an epsilon-level
// tolerance; a fixed absolute tolerance would reject valid float pairs or accept spurious double ones.
template <typename Scalar>
Scalar conjugate_tolerance(Scalar conj_tol_scale, Scalar magnitude)
{
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    return conj_tol_scale * eps * std::max(magnitude, eps);
}

template <typename Scalar, std::size_t N>
bool match_conjugate(const std::array<std::complex<Scalar>, N>& poles, std::array<bool, N>& matched, std::size_t i, Scalar conj_tol_scale)
{
    for(std::size_t j = i + 1; j < N; ++j)
    {
        if(matched[j])
            continue;
        const Scalar scale = std::max(std::abs(poles[i]), std::abs(poles[j]));
        const Scalar tol = conjugate_tolerance(conj_tol_scale, scale);
        if(std::abs(poles[i].real() - poles[j].real()) <= tol && std::abs(poles[i].imag() + poles[j].imag()) <= tol)
        {
            matched[j] = true;
            return true;
        }
    }
    return false;
}

// A real-coefficient polynomial has roots that are real or come in complex-conjugate pairs.
// @cite franklin2015 -- Franklin, Powell & Emami-Naeini, "Feedback Control of Dynamic Systems", 2015, Ch. 7
template <typename Scalar, std::size_t N>
bool validate_conjugate_pairs(const std::array<std::complex<Scalar>, N>& poles, Scalar conj_tol_scale = Scalar{2} * Scalar{N})
{
    std::array<bool, N> matched{};
    for(std::size_t i = 0; i < N; ++i)
    {
        if(matched[i] || std::abs(poles[i].imag()) <= conjugate_tolerance(conj_tol_scale, std::abs(poles[i])))
            continue;
        if(!match_conjugate(poles, matched, i, conj_tol_scale))
            return false;
    }
    return true;
}

// Returns [a_0, ..., a_{N-1}] of p(s) = s^N + a_{N-1} s^{N-1} + ... + a_0 = prod_i (s - p_i), expanded by
// repeated multiplication (Vieta's formulas); each multiplication runs from high to low degree to update in place.
// @cite strang2016 -- Strang, "Introduction to Linear Algebra", 2016, Ch. 6 (eigenvalues and characteristic polynomial)
template <typename Scalar, std::size_t N>
std::array<Scalar, N> char_poly_coeffs(const std::array<std::complex<Scalar>, N>& poles)
{
    std::array<std::complex<Scalar>, N + 1> c{};
    c[0] = std::complex<Scalar>{1, 0};

    for(std::size_t i = 0; i < N; ++i)
    {
        for(std::size_t k = i + 1; k > 0; --k)
            c[k] = c[k - 1] - poles[i] * c[k];
        c[0] = -poles[i] * c[0];
    }

    std::array<Scalar, N> result;
    for(std::size_t k = 0; k < N; ++k)
        result[k] = c[k].real();
    return result;
}

}

#endif
