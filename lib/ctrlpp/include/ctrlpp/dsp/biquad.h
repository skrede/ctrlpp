#ifndef HPP_GUARD_CTRLPP_DSP_BIQUAD_H
#define HPP_GUARD_CTRLPP_DSP_BIQUAD_H

/// @brief Second-order IIR (biquad) filter with transposed direct form II.
///
/// @cite oppenheim1997 -- Oppenheim & Willsky, "Signals and Systems", 1997
/// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae for Audio EQ Biquad Filter Coefficients", 2005

#include "ctrlpp/dsp/discrete_filter.h"

#include <array>
#include <cmath>
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

template <typename Scalar>
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
        if(std::abs(denom) < Scalar{1e-12})
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

    static auto low_pass(Scalar cutoff_hz, Scalar sample_hz) -> biquad
    {
        auto const w0 = Scalar{2} * std::numbers::pi_v<Scalar> * cutoff_hz / sample_hz;
        auto const cos_w0 = std::cos(w0);
        auto const sin_w0 = std::sin(w0);
        auto const alpha = sin_w0 / (Scalar{2} * std::numbers::sqrt2_v<Scalar>);

        auto const a0_inv = Scalar{1} / (Scalar{1} + alpha);
        return biquad{biquad_coeffs<Scalar>{
            .b0 = ((Scalar{1} - cos_w0) / Scalar{2}) * a0_inv,
            .b1 = (Scalar{1} - cos_w0) * a0_inv,
            .b2 = ((Scalar{1} - cos_w0) / Scalar{2}) * a0_inv,
            .a1 = (Scalar{-2} * cos_w0) * a0_inv,
            .a2 = (Scalar{1} - alpha) * a0_inv,
        }};
    }

    static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q) -> biquad
    {
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

    static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz) -> biquad
    {
        auto const wc = Scalar{2} * sample_hz * std::tan(std::numbers::pi_v<Scalar> * bandwidth_hz / sample_hz);
        auto const k = Scalar{2} * sample_hz;
        auto const a0_inv = Scalar{1} / (k + wc);
        return biquad{biquad_coeffs<Scalar>{
            .b0 = k * a0_inv,
            .b1 = -k * a0_inv,
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
            if(std::abs(denom) < Scalar{1e-12})
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

template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
auto make_butterworth(Scalar cutoff_hz, Scalar sample_hz) -> cascaded_biquad<Scalar, Order / 2>
{
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
/// @cite oppenheim1997 -- Oppenheim & Willsky, "Signals and Systems", 1997, Ch. 7
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
    if(std::abs(dc_gain) > Scalar{1e-15})
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

} // namespace detail

template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
auto make_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db) -> cascaded_biquad<Scalar, Order / 2>
{
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

    return cascaded_biquad<Scalar, num_sections>{sections};
}

namespace detail
{

static_assert(discrete_filter<biquad<double>>);
static_assert(discrete_filter<biquad<float>>);
static_assert(discrete_filter<cascaded_biquad<double, 2>>);

} // namespace detail

} // namespace ctrlpp

#endif
