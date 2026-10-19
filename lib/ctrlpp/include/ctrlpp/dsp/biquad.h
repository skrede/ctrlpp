#ifndef HPP_GUARD_CTRLPP_DSP_BIQUAD_H
#define HPP_GUARD_CTRLPP_DSP_BIQUAD_H

/// @brief Second-order IIR (biquad) filter with transposed direct form II.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 6 (DF-II / TDF-II structures)
/// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae for Audio EQ Biquad Filter Coefficients", 2005

#include "ctrlpp/expected.h"

#include "ctrlpp/dsp/dsp_types.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <numbers>

namespace ctrlpp
{

template <typename Scalar>
struct biquad_coeffs
{
    Scalar b0{}, b1{}, b2{};
    Scalar a1{}, a2{};
};

namespace detail
{

/// Threshold below which a biquad DC-gain quantity (the "1 + a1 + a2"
/// denominator, or an accumulated cascade DC gain) is treated as numerically
/// singular. Scaled by the instantiated Scalar's own machine epsilon and the
/// natural magnitude of the summed operands, rather than a fixed absolute
/// calibrated only for double, so the guard fires correctly at both float and
/// double precision. The default coefficient is twice the number of terms
/// entangled in the sum being tested, the same "twice the operand count"
/// convention used by the pole-placement conjugate-pair tolerance, reflecting
/// the rounding that accumulates across those terms.
template <typename Scalar>
auto biquad_singular_tol(Scalar scale, Scalar coeff = Scalar{6}) -> Scalar
{
    return coeff * std::numeric_limits<Scalar>::epsilon() * scale;
}

/// Validate the (design frequency, sample rate) pair shared by every biquad
/// design factory. Every bound is a strict domain requirement, so no epsilon
/// tolerance is involved: the sample rate of a sampled-data design must be
/// strictly positive, and the design frequency must lie strictly inside the
/// open interval (0, sample_hz / 2) because a discrete-time filter can only
/// realize a response strictly below half the sample rate; at or above it the
/// design frequency aliases (Nyquist criterion, oppenheim2010dsp Ch. 4).
/// Reports the matching `dsp_error` for a rejected pair through the library's
/// result type, and an engaged void result for a valid one.
template <typename Scalar>
auto validate_biquad_design(Scalar freq_hz, Scalar sample_hz) -> ctrlpp::expected<void, dsp_error>
{
    if(!std::isfinite(freq_hz) || !std::isfinite(sample_hz))
        return ctrlpp::unexpected(dsp_error::non_finite_input);
    if(sample_hz <= Scalar{0})
        return ctrlpp::unexpected(dsp_error::non_positive_sample_rate);
    if(freq_hz <= Scalar{0} || freq_hz >= sample_hz / Scalar{2})
        return ctrlpp::unexpected(dsp_error::cutoff_exceeds_nyquist);
    return {};
}

/// The representation error of a section's 1 + a1 + a2 with a1 stored near -2 and a2 near 1: half an ulp of a1 in
/// [1, 2) is eps / 2 and half an ulp of a2 in [1/2, 1) is eps / 4, 3/4 eps in all (oppenheim2010dsp Sec. 6.8,
/// coefficient quantization with poles near z = 1).
template <typename Scalar>
auto second_order_resolution_count() -> Scalar
{
    return Scalar{3} / Scalar{4};
}

template <typename Scalar>
struct offset_section
{
    biquad_coeffs<Scalar> coeffs;
    Scalar dc_denominator;
};

/// Unit-scale bilinear low-pass section with a1 and a2 stored as offsets from their limits -2 and 1, so the design
/// gap 1 + a1 + a2 = 4 beta0 / den0 is never formed by cancellation.
template <typename Scalar>
auto bilinear_low_pass_section(Scalar omega, Scalar g1, Scalar g0) -> offset_section<Scalar>
{
    auto const beta1 = omega * g1;
    auto const beta0 = omega * omega * g0;
    auto const den0 = Scalar{1} + beta1 + beta0;
    auto const b0 = beta0 / den0;
    return offset_section<Scalar>{
        .coeffs = biquad_coeffs<Scalar>{
            .b0 = b0,
            .b1 = Scalar{2} * b0,
            .b2 = b0,
            .a1 = Scalar{-2} + (Scalar{2} * beta1 + Scalar{4} * beta0) / den0,
            .a2 = Scalar{1} - Scalar{2} * beta1 / den0,
        },
        .dc_denominator = Scalar{4} * beta0 / den0,
    };
}

template <typename Scalar>
auto validate_resolution_factor(Scalar resolution_factor) -> ctrlpp::expected<void, dsp_error>
{
    if(!std::isfinite(resolution_factor) || !(resolution_factor > Scalar{0}))
        return ctrlpp::unexpected(dsp_error::invalid_resolution_factor);
    return {};
}

template <typename Scalar>
auto check_coefficient_resolution(Scalar dc_denominator, Scalar resolution_factor) -> ctrlpp::expected<void, dsp_error>
{
    if(!(resolution_factor * std::numeric_limits<Scalar>::epsilon() < dc_denominator))
        return ctrlpp::unexpected(dsp_error::cutoff_below_coefficient_resolution);
    return {};
}

/// Omega = tan(pi fc / fs), with the ratio formed first so pi times a large frequency cannot overflow. A result that
/// is not finite and strictly positive means the frequency rounded onto zero or onto Nyquist in Scalar arithmetic.
template <typename Scalar>
auto prewarped_frequency(Scalar freq_hz, Scalar sample_hz) -> ctrlpp::expected<Scalar, dsp_error>
{
    auto const omega = std::tan(std::numbers::pi_v<Scalar> * (freq_hz / sample_hz));
    if(!std::isfinite(omega) || !(omega > Scalar{0}))
        return ctrlpp::unexpected(dsp_error::cutoff_exceeds_nyquist);
    return omega;
}

template <typename Scalar>
auto coefficients_finite(biquad_coeffs<Scalar> const& c) -> bool
{
    return std::isfinite(c.b0) && std::isfinite(c.b1) && std::isfinite(c.b2) && std::isfinite(c.a1) && std::isfinite(c.a2);
}

template <typename Scalar>
auto validated_prewarp(Scalar freq_hz, Scalar sample_hz, Scalar resolution_factor) -> ctrlpp::expected<Scalar, dsp_error>
{
    if(auto const validated = validate_biquad_design(freq_hz, sample_hz); !validated.has_value())
        return ctrlpp::unexpected(validated.error());
    if(auto const factor = validate_resolution_factor(resolution_factor); !factor.has_value())
        return ctrlpp::unexpected(factor.error());
    return prewarped_frequency(freq_hz, sample_hz);
}

/// A section whose design gap is not finite cannot be judged against the floor, so finiteness is refused first.
template <typename Scalar>
auto resolved_coefficients(offset_section<Scalar> const& section, Scalar resolution_factor) -> ctrlpp::expected<biquad_coeffs<Scalar>, dsp_error>
{
    if(!coefficients_finite(section.coeffs) || !std::isfinite(section.dc_denominator))
        return ctrlpp::unexpected(dsp_error::non_finite_input);
    if(auto const resolved = check_coefficient_resolution(section.dc_denominator, resolution_factor); !resolved.has_value())
        return ctrlpp::unexpected(resolved.error());
    return section.coeffs;
}

}

template <ctrlpp_floating_scalar Scalar>
class biquad
{
public:
    using scalar_type = Scalar;

    constexpr biquad() = default;
    explicit constexpr biquad(biquad_coeffs<Scalar> c) : c_{c} {}

    auto process(Scalar x) -> Scalar
    {
        auto const y = c_.b0 * x + w_[0];
        w_[0] = c_.b1 * x - c_.a1 * y + w_[1];
        w_[1] = c_.b2 * x - c_.a2 * y;
        return y;
    }

    void reset() { w_.fill(Scalar{0}); }

    void reset(Scalar value)
    {
        auto const denom = Scalar{1} + c_.a1 + c_.a2;
        auto const denom_scale = Scalar{1} + std::abs(c_.a1) + std::abs(c_.a2);
        if(std::abs(denom) < detail::biquad_singular_tol(denom_scale))
        {
            w_.fill(Scalar{0});
            return;
        }
        auto const dc_gain = (c_.b0 + c_.b1 + c_.b2) / denom;
        auto const y_ss = value * dc_gain;
        w_[0] = value * (c_.b1 + c_.b2) - y_ss * (c_.a1 + c_.a2);
        w_[1] = value * c_.b2 - y_ss * c_.a2;
    }

    auto coefficients() const -> biquad_coeffs<Scalar> const& { return c_; }

    /// Second-order Butterworth low-pass biquad section via the RBJ cookbook
    /// formulas (analog-prototype design mapped through the bilinear transform),
    /// formed at unit scale in Omega = tan(pi cutoff_hz / sample_hz).
    ///
    /// Rejections, checked in order, each with the matching `dsp_error`:
    /// non-finite parameters; sample_hz <= 0; cutoff_hz outside the open
    /// interval (0, sample_hz / 2) (Nyquist criterion, oppenheim2010dsp Ch. 4);
    /// a resolution_factor that is not finite and strictly positive; a cutoff
    /// ratio that rounds onto zero or Nyquist in the pre-warp
    /// (cutoff_exceeds_nyquist); a non-finite coefficient; a design gap
    /// 1 + a1 + a2 at or below resolution_factor * eps
    /// (cutoff_below_coefficient_resolution).
    ///
    /// resolution_factor defaults to the counted representation error of
    /// 1 + a1 + a2, so the default refuses a design whose gap the coefficients
    /// cannot carry to one correct digit; a caller wanting relative accuracy tau
    /// passes that count divided by tau.
    ///
    /// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (LPF section)
    /// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (bilinear transform of analog prototypes)
    static auto low_pass(Scalar cutoff_hz, Scalar sample_hz, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
        -> expected<biquad, dsp_error>
    {
        auto const omega = detail::validated_prewarp(cutoff_hz, sample_hz, resolution_factor);
        if(!omega.has_value())
            return unexpected(omega.error());
        // RBJ's Butterworth low-pass rewritten exactly in Omega: Q = 1 / sqrt(2) is g1 = sqrt(2).
        auto const section = detail::resolved_coefficients(detail::bilinear_low_pass_section(*omega, std::numbers::sqrt2_v<Scalar>, Scalar{1}), resolution_factor);
        if(!section.has_value())
            return unexpected(section.error());
        return biquad{*section};
    }

    /// Notch (band-stop) biquad section via the RBJ cookbook formulas.
    ///
    /// Rejects non-finite parameters, sample_hz <= 0, freq_hz outside the open
    /// interval (0, sample_hz / 2) (Nyquist criterion, oppenheim2010dsp Ch. 4),
    /// and q <= 0 (the design bandwidth alpha = sin(w0) / (2 q) requires a
    /// strictly positive quality factor) with the matching `dsp_error`.
    ///
    /// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (notch section)
    /// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 6
    static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q) -> expected<biquad, dsp_error>
    {
        if(!std::isfinite(q))
            return unexpected(dsp_error::non_finite_input);
        if(auto const validated = detail::validate_biquad_design(freq_hz, sample_hz); !validated.has_value())
            return unexpected(validated.error());
        if(q <= Scalar{0})
            return unexpected(dsp_error::non_positive_q);

        auto const w0 = Scalar{2} * std::numbers::pi_v<Scalar> * freq_hz / sample_hz;
        auto const cos_w0 = std::cos(w0);
        auto const alpha = std::sin(w0) / (Scalar{2} * q);

        auto const a0_inv = Scalar{1} / (Scalar{1} + alpha);
        return biquad{biquad_coeffs<Scalar>{
            .b0 = a0_inv,
            .b1 = (Scalar{-2} * cos_w0) * a0_inv,
            .b2 = a0_inv,
            .a1 = (Scalar{-2} * cos_w0) * a0_inv,
            .a2 = (Scalar{1} - alpha) * a0_inv,
        }};
    }

    /// "Dirty" derivative: a band-limited differentiator with analog prototype
    /// H(s) = wc*s / (s + wc), discretized by the bilinear transform with
    /// frequency pre-warping. The wc numerator factor makes this a true
    /// differentiator (|H| -> 2*pi*f) up to the cutoff wc, above which it rolls
    /// off; without it the section would be a unity-gain high-pass whose
    /// passband gain is wrong by 1/wc. The bilinear map with k = 2*fs gives
    /// numerator coefficients b0 = wc*k/(k+wc), b1 = -wc*k/(k+wc); only the
    /// numerator carries the differentiator gain, so the denominator is
    /// unchanged.
    ///
    /// Rejects non-finite parameters, sample_hz <= 0, and bandwidth_hz outside
    /// the open interval (0, sample_hz / 2) (Nyquist criterion, oppenheim2010dsp
    /// Ch. 4; at half the sample rate the pre-warped cutoff tan(pi/2) also
    /// diverges) with the matching `dsp_error`.
    ///
    /// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (bilinear transform with pre-warping)
    static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz) -> expected<biquad, dsp_error>
    {
        if(auto const validated = detail::validate_biquad_design(bandwidth_hz, sample_hz); !validated.has_value())
            return unexpected(validated.error());

        auto const wc = Scalar{2} * sample_hz * std::tan(std::numbers::pi_v<Scalar> * bandwidth_hz / sample_hz);
        auto const k = Scalar{2} * sample_hz;
        auto const a0_inv = Scalar{1} / (k + wc);
        return biquad{biquad_coeffs<Scalar>{
            .b0 = wc * k * a0_inv,
            .b1 = -wc * k * a0_inv,
            .b2 = Scalar{0},
            .a1 = (wc - k) * a0_inv,
            .a2 = Scalar{0},
        }};
    }

private:
    biquad_coeffs<Scalar> c_{};
    std::array<Scalar, 2> w_{};
};

template <typename Scalar>
biquad(biquad_coeffs<Scalar>) -> biquad<Scalar>;

template <typename Scalar, std::size_t N>
    requires(N >= 1)
class cascaded_biquad
{
public:
    using scalar_type = Scalar;

    explicit constexpr cascaded_biquad(std::array<biquad<Scalar>, N> sections) : sections_{sections} {}

    auto process(Scalar x) -> Scalar
    {
        for(auto& s : sections_)
        {
            x = s.process(x);
        }
        return x;
    }

    void reset()
    {
        for(auto& s : sections_)
        {
            s.reset();
        }
    }

    void reset(Scalar value)
    {
        for(auto& s : sections_)
        {
            s.reset(value);
            auto const& c = s.coefficients();
            auto const denom = Scalar{1} + c.a1 + c.a2;
            auto const denom_scale = Scalar{1} + std::abs(c.a1) + std::abs(c.a2);
            if(std::abs(denom) < detail::biquad_singular_tol(denom_scale))
            {
                value = Scalar{0};
            }
            else
            {
                value *= (c.b0 + c.b1 + c.b2) / denom;
            }
        }
    }

    auto section(std::size_t i) const -> biquad<Scalar> const& { return sections_[i]; }

    auto section(std::size_t i) -> biquad<Scalar>& { return sections_[i]; }

private:
    std::array<biquad<Scalar>, N> sections_{};
};

/// Cascade of second-order Butterworth low-pass sections.
///
/// Builds an even-order Butterworth filter as a product of biquad sections
/// from the analog-prototype pole locations on the unit circle, mapped through
/// the bilinear transform with frequency pre-warping. Section k is the unit-scale
/// section with g1 = 2 cos(theta_k) = 1 / q_k.
///
/// Rejects, in the order and with the errors `biquad::low_pass` documents,
/// every section's design gap checked against the coefficient-resolution floor
/// (the smallest gap decides).
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (analog-prototype Butterworth, bilinear pre-warp)
/// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (per-section LPF coefficients)
template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
auto make_butterworth(Scalar cutoff_hz, Scalar sample_hz, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
    -> expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>
{
    auto const omega = detail::validated_prewarp(cutoff_hz, sample_hz, resolution_factor);
    if(!omega.has_value())
        return unexpected(omega.error());

    std::array<biquad<Scalar>, Order / 2> sections{};
    for(std::size_t k = 0; k < Order / 2; ++k)
    {
        auto const theta = std::numbers::pi_v<Scalar> * static_cast<Scalar>(2 * k + 1) / static_cast<Scalar>(2 * Order);
        auto const section = detail::resolved_coefficients(detail::bilinear_low_pass_section(*omega, Scalar{2} * std::cos(theta), Scalar{1}), resolution_factor);
        if(!section.has_value())
            return unexpected(section.error());
        sections[k] = biquad{*section};
    }
    return cascaded_biquad<Scalar, Order / 2>{sections};
}

namespace detail
{

/// @brief One Chebyshev Type I section: the unit-scale bilinear section of the prototype pole pair sigma +- j imag,
/// g1 = -2 sigma and g0 = sigma^2 + imag^2, which is the bilinear map divided through by (2 fs)^2.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I analog prototype, bilinear pre-warp)
template <typename Scalar, std::size_t Order>
auto chebyshev1_section(std::size_t k, Scalar sinh_v, Scalar cosh_v, Scalar omega, Scalar resolution_factor)
    -> ctrlpp::expected<biquad_coeffs<Scalar>, dsp_error>
{
    auto const theta = std::numbers::pi_v<Scalar> * static_cast<Scalar>(2 * k + 1) / static_cast<Scalar>(2 * Order);
    auto const sigma = -sinh_v * std::sin(theta);
    auto const imag = cosh_v * std::cos(theta);
    return resolved_coefficients(bilinear_low_pass_section(omega, Scalar{-2} * sigma, sigma * sigma + imag * imag), resolution_factor);
}

/// @brief Normalize cascade DC gain so max passband gain = 0 dB.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I passband normalization)
/// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005
template <typename Scalar, std::size_t N>
void normalize_chebyshev1_dc(std::array<biquad<Scalar>, N>& sections, Scalar eps)
{
    auto dc_gain = Scalar{1};
    for(auto const& s : sections)
    {
        auto const& c = s.coefficients();
        dc_gain *= (c.b0 + c.b1 + c.b2) / (Scalar{1} + c.a1 + c.a2);
    }
    auto const target_dc = Scalar{1} / std::sqrt(Scalar{1} + eps * eps);
    if(std::abs(dc_gain) > biquad_singular_tol(target_dc, Scalar{2} * Scalar{N}))
    {
        auto const correction = target_dc / dc_gain;
        auto const old = sections[0].coefficients();
        sections[0] = biquad{biquad_coeffs<Scalar>{
            .b0 = old.b0 * correction,
            .b1 = old.b1 * correction,
            .b2 = old.b2 * correction,
            .a1 = old.a1,
            .a2 = old.a2,
        }};
    }
}

/// True when every coefficient of every section of a designed cascade is
/// finite. A design factory sweeps its emitted sections with this before
/// reporting success, so a successful design never hands back a filter whose
/// difference equation immediately contaminates its state with NaN or infinity.
/// The predicate is exact; there is no threshold involved.
template <typename Scalar, std::size_t N>
auto sections_all_finite(std::array<biquad<Scalar>, N> const& sections) -> bool
{
    for(auto const& s : sections)
    {
        if(!coefficients_finite(s.coefficients()))
            return false;
    }
    return true;
}

}

/// Cascade of second-order Chebyshev Type I low-pass sections.
///
/// Builds an even-order Chebyshev Type I filter as a product of biquad
/// sections from the analog-prototype poles on an ellipse parameterized by
/// the passband ripple, mapped through the bilinear transform with frequency
/// pre-warping and formed at unit scale. The DC gain is normalized so the max
/// passband gain is 0 dB.
///
/// Rejections, checked in order:
///  * NaN or infinite ripple_db -> dsp_error::non_finite_input
///  * ripple_db <= 0            -> dsp_error::non_positive_ripple. The ripple
///    factor is eps = sqrt(10^(ripple_db / 10) - 1); its radicand is
///    non-positive for every ripple_db <= 0, and at exactly zero the following
///    asinh(1 / eps) takes an infinite argument. That is the exact domain of
///    the design, not a tolerance, so no constant enters the bound.
///  * non-finite cutoff_hz or sample_hz -> dsp_error::non_finite_input
///  * sample_hz <= 0            -> dsp_error::non_positive_sample_rate
///  * cutoff_hz outside the open interval (0, sample_hz / 2)
///    -> dsp_error::cutoff_exceeds_nyquist (Nyquist criterion,
///    oppenheim2010dsp Ch. 4)
///  * resolution_factor not finite and strictly positive
///    -> dsp_error::invalid_resolution_factor
///  * a cutoff ratio that rounds onto zero or Nyquist in the pre-warp
///    -> dsp_error::cutoff_exceeds_nyquist
///  * per section, a non-finite coefficient -> dsp_error::non_finite_input,
///    then a design gap 1 + a1 + a2 at or below resolution_factor * eps
///    -> dsp_error::cutoff_below_coefficient_resolution (the smallest gap
///    decides)
///  * a non-finite emitted coefficient after the DC normalization
///    -> dsp_error::non_finite_input. The emitted sections are swept before
///    the success return, so a successful design never carries a coefficient
///    a filter cannot run.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I IIR design)
template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
auto make_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
    -> expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>
{
    if(!std::isfinite(ripple_db))
        return unexpected(dsp_error::non_finite_input);
    if(ripple_db <= Scalar{0})
        return unexpected(dsp_error::non_positive_ripple);
    auto const omega = detail::validated_prewarp(cutoff_hz, sample_hz, resolution_factor);
    if(!omega.has_value())
        return unexpected(omega.error());

    // 10^(ripple_db / 10) - 1 through expm1, so a small ripple does not cancel against the 1.
    auto const eps = std::sqrt(std::expm1(std::numbers::ln10_v<Scalar> * ripple_db / Scalar{10}));
    auto const v = std::asinh(Scalar{1} / eps) / static_cast<Scalar>(Order);

    std::array<biquad<Scalar>, Order / 2> sections{};
    for(std::size_t k = 0; k < Order / 2; ++k)
    {
        auto const section = detail::chebyshev1_section<Scalar, Order>(k, std::sinh(v), std::cosh(v), *omega, resolution_factor);
        if(!section.has_value())
            return unexpected(section.error());
        sections[k] = biquad{*section};
    }
    detail::normalize_chebyshev1_dc(sections, eps);
    if(!detail::sections_all_finite(sections))
        return unexpected(dsp_error::non_finite_input);
    return cascaded_biquad<Scalar, Order / 2>{sections};
}

namespace detail
{

static_assert(discrete_filter<biquad<double>>);
static_assert(discrete_filter<biquad<float>>);
static_assert(discrete_filter<cascaded_biquad<double, 2>>);

}

}

#endif
