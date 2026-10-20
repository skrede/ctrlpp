#ifndef HPP_GUARD_CTRLPP_MODEL_CONVERSION_H
#define HPP_GUARD_CTRLPP_MODEL_CONVERSION_H

/// @brief Transfer-function to/from state-space conversion (controllable canonical form; characteristic
/// polynomials from the eigenvalues of a diagonally balanced matrix).
///
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Ch. 2-3
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.5 (QR algorithm, backward stable)
/// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 algorithm (balancing)

#include "ctrlpp/expected.h"

#include "ctrlpp/model/state_space.h"
#include "ctrlpp/model/conversion_types.h"
#include "ctrlpp/model/transfer_function.h"

#include "ctrlpp/detail/characteristic_polynomial.h"

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t N>
auto coefficients_finite(const std::array<Scalar, N>& coefficients) -> bool
{
    return std::ranges::all_of(coefficients, [](Scalar c) { return std::isfinite(c); });
}

template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
auto check_transfer_function(const transfer_function<Scalar, NumDeg, DenDeg>& tf) -> ctrlpp::expected<void, conversion_error>
{
    if(!coefficients_finite(tf.numerator) || !coefficients_finite(tf.denominator))
        return ctrlpp::unexpected(conversion_error::non_finite_input);
    if(tf.denominator[0] == Scalar{0})
        return ctrlpp::unexpected(conversion_error::zero_leading_coefficient);
    return {};
}

template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
auto monic_padded(const transfer_function<Scalar, NumDeg, DenDeg>& tf) -> transfer_function<Scalar, DenDeg, DenDeg>
{
    const Scalar a0 = tf.denominator[0];
    transfer_function<Scalar, DenDeg, DenDeg> monic{};
    for(std::size_t i = 0; i <= DenDeg; ++i)
        monic.denominator[i] = tf.denominator[i] / a0;
    for(std::size_t i = 0; i <= NumDeg; ++i)
        monic.numerator[DenDeg - NumDeg + i] = tf.numerator[i] / a0;
    return monic;
}

template <bool Biproper, typename Scalar, std::size_t N>
auto companion_realization(const transfer_function<Scalar, N, N>& monic) -> continuous_state_space<Scalar, N, 1, 1>
{
    constexpr int n = static_cast<int>(N);
    continuous_state_space<Scalar, N, 1, 1> sys{};
    for(int i = 0; i + 1 < n; ++i)
        sys.A(i, i + 1) = Scalar{1};
    sys.B(n - 1, 0) = Scalar{1};
    for(int j = 0; j < n; ++j)
    {
        const std::size_t k = N - static_cast<std::size_t>(j);
        sys.A(n - 1, j) = -monic.denominator[k];
        if constexpr(Biproper)
            sys.C(0, j) = monic.numerator[k] - monic.numerator[0] * monic.denominator[k];
        else
            sys.C(0, j) = monic.numerator[k];
    }
    if constexpr(Biproper)
        sys.D(0, 0) = monic.numerator[0];
    return sys;
}

/// 2^round(log2(||A|| / (||B|| ||C||))), formed from the logarithms so no product of norms overflows and
/// clamped to the normal exponent range so alpha stays a finite, exact power of two.
template <typename Scalar>
auto rank_one_scale(Scalar norm_a, Scalar norm_b, Scalar norm_c) -> Scalar
{
    if(norm_a == Scalar{0})
        return Scalar{1};
    constexpr int lowest = std::numeric_limits<Scalar>::min_exponent - 1;
    constexpr int highest = std::numeric_limits<Scalar>::max_exponent - 1;
    const auto e = static_cast<int>(std::round(std::log2(norm_a) - std::log2(norm_b) - std::log2(norm_c)));
    return std::ldexp(Scalar{1}, std::clamp(e, lowest, highest));
}

/// The rank-one identity det(sI - A + alpha B C) = det(sI - A) (1 + alpha C (sI - A)^{-1} B) makes the
/// strictly proper numerator (poly(A - alpha B C) - poly(A)) / alpha for every alpha != 0.
template <typename Scalar, std::size_t NX>
auto transfer_numerator(const continuous_state_space<Scalar, NX, 1, 1>& sys, const std::array<Scalar, NX + 1>& den)
    -> ctrlpp::expected<std::array<Scalar, NX + 1>, conversion_error>
{
    const Scalar d = sys.D(0, 0);
    std::array<Scalar, NX + 1> num{};
    std::ranges::transform(den, num.begin(), [d](Scalar c) { return d * c; });
    const Scalar norm_b = sys.B.stableNorm();
    const Scalar norm_c = sys.C.stableNorm();
    if(norm_b == Scalar{0} || norm_c == Scalar{0})
        return num;
    // Eigen 3.4.0's stableNorm asserts on a fixed-size non-vector matrix; the reshaped view takes the vector path.
    const Scalar alpha = rank_one_scale(sys.A.reshaped().stableNorm(), norm_b, norm_c);
    const auto shifted = characteristic_polynomial<Scalar, int(NX)>(sys.A - alpha * sys.B * sys.C);
    if(!shifted)
        return ctrlpp::unexpected(shifted.error());
    for(std::size_t k = 1; k <= NX; ++k)
        num[k] = ((*shifted)[k] - den[k]) / alpha + d * den[k];
    return num;
}

}

/// @brief Transfer function to controllable canonical form state-space.
///
/// H(s) = num(s) / den(s), coefficients highest-degree-first, normalized to a
/// monic denominator. Returns continuous_state_space with NX = DenDeg states.
///
/// Refusals, checked in order:
///  * a NaN or infinite coefficient -> conversion_error::non_finite_input.
///  * denominator[0] == 0 -> conversion_error::zero_leading_coefficient.
///  * a realization entry overflows (a tiny leading coefficient's quotient)
///    -> conversion_error::unrepresentable_result.
///
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.1 (controllable canonical form)
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015, Ch. 7
template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
    requires(NumDeg <= DenDeg)
auto tf2ss(const transfer_function<Scalar, NumDeg, DenDeg>& tf)
    -> ctrlpp::expected<continuous_state_space<Scalar, DenDeg, 1, 1>, conversion_error>
{
    static_assert(DenDeg >= 1, "Denominator degree must be at least 1");
    if(auto const gate = detail::check_transfer_function(tf); !gate)
        return ctrlpp::unexpected(gate.error());
    auto const sys = detail::companion_realization<NumDeg == DenDeg>(detail::monic_padded(tf));
    if(!sys.A.allFinite() || !sys.C.allFinite() || !sys.D.allFinite())
        return ctrlpp::unexpected(conversion_error::unrepresentable_result);
    return sys;
}

/// @brief State-space to transfer function, H(s) = C (sI - A)^{-1} B + D, for SISO systems.
///
/// Returns transfer_function<Scalar, NX, NX>, coefficients highest-degree-first. The denominator is the
/// characteristic polynomial of A from the eigenvalues of a diagonally balanced copy; the numerator is
/// D poly(A) plus (poly(A - alpha B C) - poly(A)) / alpha through the same routine, with alpha the power
/// of two nearest ||A||_F / (||B|| ||C||). The leading numerator coefficient is D exactly; B = 0 or C = 0
/// gives D poly(A).
///
/// Refusals, checked in order:
///  * a NaN or infinite entry in A, B, C or D -> conversion_error::non_finite_input.
///  * the eigensolver does not converge -> conversion_error::eigenvalues_not_converged.
///  * a coefficient is not finite -> conversion_error::unrepresentable_result.
template <typename Scalar, std::size_t NX>
auto ss2tf(const continuous_state_space<Scalar, NX, 1, 1>& sys)
    -> ctrlpp::expected<transfer_function<Scalar, NX, NX>, conversion_error>
{
    if(!sys.A.allFinite() || !sys.B.allFinite() || !sys.C.allFinite() || !sys.D.allFinite())
        return ctrlpp::unexpected(conversion_error::non_finite_input);
    const auto den = detail::characteristic_polynomial<Scalar, int(NX)>(sys.A);
    if(!den)
        return ctrlpp::unexpected(den.error());
    const auto num = detail::transfer_numerator(sys, *den);
    if(!num)
        return ctrlpp::unexpected(num.error());
    if(!detail::coefficients_finite(*num) || !detail::coefficients_finite(*den))
        return ctrlpp::unexpected(conversion_error::unrepresentable_result);
    return transfer_function<Scalar, NX, NX>{*num, *den};
}

}

#endif
