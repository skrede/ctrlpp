#ifndef HPP_GUARD_CTRLPP_DSP_CASCADED_BIQUAD_H
#define HPP_GUARD_CTRLPP_DSP_CASCADED_BIQUAD_H

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/dsp_types.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/dsp/detail/biquad_design.h"
#include "ctrlpp/dsp/detail/biquad_sections.h"

#include "ctrlpp/expected.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace ctrlpp
{

template <typename Scalar, std::size_t N>
    requires(N >= 1)
class cascaded_biquad
{
public:
    using scalar_type = Scalar;

    explicit constexpr cascaded_biquad(std::array<biquad<Scalar>, N> sections)
        : m_sections{sections}
    {
    }

    Scalar process(Scalar x)
    {
        for(auto& s : m_sections)
        {
            x = s.process(x);
        }
        return x;
    }

    void reset()
    {
        for(auto& s : m_sections)
        {
            s.reset();
        }
    }

    void reset(Scalar value)
    {
        for(auto& s : m_sections)
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

    biquad<Scalar> const& section(std::size_t i) const { return m_sections[i]; }

    biquad<Scalar>& section(std::size_t i) { return m_sections[i]; }

private:
    std::array<biquad<Scalar>, N> m_sections;
};

// Section k is the unit-scale section with g1 = 2 cos(theta_k) = 1 / q_k.
// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (analog-prototype Butterworth, bilinear pre-warp)
// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (per-section LPF coefficients)
template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
expected<cascaded_biquad<Scalar, Order / 2>, dsp_error> make_butterworth(Scalar cutoff_hz, Scalar sample_hz, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
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

// An even-order Chebyshev Type I response starts at its ripple trough, so the DC gain is set to 1 / sqrt(1 + eps^2)
// for a 0 dB passband peak.
// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I passband normalization)
// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005
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
        sections[0] = biquad{biquad_coeffs<Scalar>{old.b0 * correction, old.b1 * correction, old.b2 * correction, old.a1, old.a2}};
    }
}

template <typename Scalar, std::size_t N>
bool sections_all_finite(std::array<biquad<Scalar>, N> const& sections)
{
    for(auto const& s : sections)
    {
        if(!coefficients_finite(s.coefficients()))
            return false;
    }
    return true;
}

}

// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I IIR design)
template <std::size_t Order, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
expected<cascaded_biquad<Scalar, Order / 2>, dsp_error> make_chebyshev1(
        Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
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

static_assert(discrete_filter<cascaded_biquad<double, 2>>);

}

}

#endif
