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
#include <optional>

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
/// Returns the matching `dsp_error` for a rejected pair, or an empty optional
/// for a valid one.
template <typename Scalar>
auto validate_biquad_design(Scalar freq_hz, Scalar sample_hz) -> std::optional<dsp_error>
{
    if(!std::isfinite(freq_hz) || !std::isfinite(sample_hz))
        return dsp_error::non_finite_input;
    if(sample_hz <= Scalar{0})
        return dsp_error::non_positive_sample_rate;
    if(freq_hz <= Scalar{0} || freq_hz >= sample_hz / Scalar{2})
        return dsp_error::cutoff_exceeds_nyquist;
    return std::nullopt;
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

    [[nodiscard]] auto coefficients() const -> biquad_coeffs<Scalar> const& { return c_; }

    /// Second-order Butterworth low-pass biquad section via the RBJ cookbook
    /// formulas (analog-prototype design mapped through the bilinear transform).
    ///
    /// Rejects non-finite parameters, sample_hz <= 0, and cutoff_hz outside the
    /// open interval (0, sample_hz / 2) (Nyquist criterion, oppenheim2010dsp
    /// Ch. 4) with the matching `dsp_error`.
    ///
    /// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (LPF section)
    /// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (bilinear transform of analog prototypes)
    [[nodiscard]] static auto low_pass(Scalar cutoff_hz, Scalar sample_hz) -> expected<biquad, dsp_error>
    {
        if(auto const err = detail::validate_biquad_design(cutoff_hz, sample_hz))
            return unexpected(*err);

        auto const w0 = Scalar{2} * std::numbers::pi_v<Scalar> * cutoff_hz / sample_hz;
        auto const cos_w0 = std::cos(w0);
        auto const sin_w0 = std::sin(w0);
        // Butterworth (maximally flat) response fixes the quality factor at
        // Q = 1/sqrt(2), so the RBJ cookbook width alpha = sin(w0)/(2*Q)
        // reduces to sin(w0)/sqrt(2). This matches make_butterworth<2> in this
        // same file (q_k = 1/sqrt(2)) and yields -3.01 dB at the cutoff with no
        // passband peaking.
        auto const alpha = sin_w0 / std::numbers::sqrt2_v<Scalar>;

        auto const a0_inv = Scalar{1} / (Scalar{1} + alpha);
        return biquad{biquad_coeffs<Scalar>{
            .b0 = ((Scalar{1} - cos_w0) / Scalar{2}) * a0_inv,
            .b1 = (Scalar{1} - cos_w0) * a0_inv,
            .b2 = ((Scalar{1} - cos_w0) / Scalar{2}) * a0_inv,
            .a1 = (Scalar{-2} * cos_w0) * a0_inv,
            .a2 = (Scalar{1} - alpha) * a0_inv,
        }};
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
    [[nodiscard]] static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q) -> expected<biquad, dsp_error>
    {
        if(!std::isfinite(q))
            return unexpected(dsp_error::non_finite_input);
        if(auto const err = detail::validate_biquad_design(freq_hz, sample_hz))
            return unexpected(*err);
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
    /// H(s) = wc*s / (s + wc), discretised by the bilinear transform with
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
    [[nodiscard]] static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz) -> expected<biquad, dsp_error>
    {
        if(auto const err = detail::validate_biquad_design(bandwidth_hz, sample_hz))
            return unexpected(*err);

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

    [[nodiscard]] auto section(std::size_t i) const -> biquad<Scalar> const& { return sections_[i]; }

    auto section(std::size_t i) -> biquad<Scalar>& { return sections_[i]; }

private:
    std::array<biquad<Scalar>, N> sections_{};
};

/// Cascade of second-order Butterworth low-pass sections.
///
/// Builds an even-order Butterworth filter as a product of biquad sections
/// from the analog-prototype pole locations on the unit circle, mapped through
/// the bilinear transform with frequency pre-warping.
///
/// Rejects non-finite parameters, sample_hz <= 0, and cutoff_hz outside the
/// open interval (0, sample_hz / 2) (Nyquist criterion, oppenheim2010dsp
/// Ch. 4) with the matching `dsp_error`.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (analog-prototype Butterworth, bilinear pre-warp)
/// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (per-section LPF coefficients)
template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
[[nodiscard]] auto make_butterworth(Scalar cutoff_hz, Scalar sample_hz)
    -> expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>
{
    if(auto const err = detail::validate_biquad_design(cutoff_hz, sample_hz))
        return unexpected(*err);

    constexpr auto num_sections = Order / 2;
    auto const w0 = Scalar{2} * std::numbers::pi_v<Scalar> * cutoff_hz / sample_hz;
    auto const cos_w0 = std::cos(w0);
    auto const sin_w0 = std::sin(w0);

    std::array<biquad<Scalar>, num_sections> sections{};
    for(std::size_t k = 0; k < num_sections; ++k)
    {
        auto const theta = std::numbers::pi_v<Scalar> * static_cast<Scalar>(2 * k + 1) / static_cast<Scalar>(2 * Order);
        auto const q_k = Scalar{1} / (Scalar{2} * std::cos(theta));
        auto const alpha = sin_w0 / (Scalar{2} * q_k);

        auto const a0_inv = Scalar{1} / (Scalar{1} + alpha);
        sections[k] = biquad{biquad_coeffs<Scalar>{
            .b0 = ((Scalar{1} - cos_w0) / Scalar{2}) * a0_inv,
            .b1 = (Scalar{1} - cos_w0) * a0_inv,
            .b2 = ((Scalar{1} - cos_w0) / Scalar{2}) * a0_inv,
            .a1 = (Scalar{-2} * cos_w0) * a0_inv,
            .a2 = (Scalar{1} - alpha) * a0_inv,
        }};
    }
    return cascaded_biquad<Scalar, num_sections>{sections};
}

namespace detail
{

/// @brief Compute one Chebyshev Type I biquad section via bilinear transform.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I analog prototype, bilinear pre-warp)
template <typename Scalar, std::size_t Order>
auto chebyshev1_section(std::size_t k, Scalar sinh_v, Scalar cosh_v, Scalar wc, Scalar sample_hz) -> biquad<Scalar>
{
    auto const theta = std::numbers::pi_v<Scalar> * static_cast<Scalar>(2 * k + 1) / static_cast<Scalar>(2 * Order);
    auto const sigma = -sinh_v * std::sin(theta);
    auto const omega = cosh_v * std::cos(theta);

    auto const pole_mag_sq = sigma * sigma + omega * omega;
    auto const b1_a = Scalar{-2} * sigma * wc;
    auto const b0_a = pole_mag_sq * wc * wc;

    auto const k_bt = Scalar{2} * sample_hz;
    auto const k2 = k_bt * k_bt;

    auto const den_z0 = k2 + b1_a * k_bt + b0_a;
    auto const a0_inv = Scalar{1} / den_z0;

    return biquad{biquad_coeffs<Scalar>{
        .b0 = b0_a * a0_inv,
        .b1 = Scalar{2} * b0_a * a0_inv,
        .b2 = b0_a * a0_inv,
        .a1 = (Scalar{2} * b0_a - Scalar{2} * k2) * a0_inv,
        .a2 = (k2 - b1_a * k_bt + b0_a) * a0_inv,
    }};
}

/// @brief Normalize cascade DC gain so max passband gain = 0 dB.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I passband normalisation)
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
        auto const& c = s.coefficients();
        if(!std::isfinite(c.b0) || !std::isfinite(c.b1) || !std::isfinite(c.b2)
           || !std::isfinite(c.a1) || !std::isfinite(c.a2))
        {
            return false;
        }
    }
    return true;
}

}

/// Cascade of second-order Chebyshev Type I low-pass sections.
///
/// Builds an even-order Chebyshev Type I filter as a product of biquad
/// sections from the analog-prototype poles on an ellipse parameterised by
/// the passband ripple, mapped through the bilinear transform with frequency
/// pre-warping. The DC gain is normalised so the max passband gain is 0 dB.
///
/// Rejections, checked in order:
///  * NaN or infinite ripple_db -> dsp_error::non_finite_input
///  * ripple_db <= 0            -> dsp_error::non_positive_ripple. The ripple
///    factor is eps = sqrt(10^(ripple_db / 10) - 1); its radicand is
///    non-positive for every ripple_db <= 0, and at exactly zero the following
///    asinh(1 / eps) takes an infinite argument. That is the exact domain of
///    the design, not a tolerance, so no constant enters the bound.
///  * sample_hz <= 0            -> dsp_error::non_positive_sample_rate
///  * cutoff_hz outside the open interval (0, sample_hz / 2)
///    -> dsp_error::cutoff_exceeds_nyquist (Nyquist criterion,
///    oppenheim2010dsp Ch. 4)
///  * a non-finite emitted coefficient -> dsp_error::non_finite_input. The
///    emitted sections are swept before the success return, so a successful
///    design never carries a coefficient a filter cannot run.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I IIR design)
template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
[[nodiscard]] auto make_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db)
    -> expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>
{
    if(!std::isfinite(ripple_db))
        return unexpected(dsp_error::non_finite_input);
    if(ripple_db <= Scalar{0})
        return unexpected(dsp_error::non_positive_ripple);
    if(auto const err = detail::validate_biquad_design(cutoff_hz, sample_hz))
        return unexpected(*err);

    constexpr auto num_sections = Order / 2;

    auto const eps = std::sqrt(std::pow(Scalar{10}, ripple_db / Scalar{10}) - Scalar{1});
    auto const v = std::asinh(Scalar{1} / eps) / static_cast<Scalar>(Order);
    auto const sinh_v = std::sinh(v);
    auto const cosh_v = std::cosh(v);
    auto const wc = Scalar{2} * sample_hz * std::tan(std::numbers::pi_v<Scalar> * cutoff_hz / sample_hz);

    std::array<biquad<Scalar>, num_sections> sections{};
    for(std::size_t k = 0; k < num_sections; ++k)
        sections[k] = detail::chebyshev1_section<Scalar, Order>(k, sinh_v, cosh_v, wc, sample_hz);

    detail::normalize_chebyshev1_dc(sections, eps);

    if(!detail::sections_all_finite(sections))
        return unexpected(dsp_error::non_finite_input);

    return cascaded_biquad<Scalar, num_sections>{sections};
}

namespace detail
{

static_assert(discrete_filter<biquad<double>>);
static_assert(discrete_filter<biquad<float>>);
static_assert(discrete_filter<cascaded_biquad<double, 2>>);

}

}

#endif
