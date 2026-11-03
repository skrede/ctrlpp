#ifndef HPP_GUARD_CTRLPP_DSP_DETAIL_BIQUAD_SECTIONS_H
#define HPP_GUARD_CTRLPP_DSP_DETAIL_BIQUAD_SECTIONS_H

// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae for Audio EQ Biquad Filter Coefficients", 2005

#include "ctrlpp/dsp/detail/biquad_design.h"

#include "ctrlpp/dsp/dsp_types.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <cstddef>
#include <numbers>

namespace ctrlpp::detail
{

// pole_offset is 1 minus the product of the poles (1 - a2, or 1 + a1 for a first-order section), kept so an
// offset-form numerator reuses it instead of re-forming it.
template <typename Scalar>
struct offset_section
{
    biquad_coeffs<Scalar> coeffs;
    Scalar dc_denominator;
    Scalar pole_offset;
};

// Unit-scale bilinear low-pass section with a1 and a2 stored as offsets from their limits -2 and 1, so the design
// gap 1 + a1 + a2 = 4 beta0 / den0 is never formed by cancellation.
template <typename Scalar>
offset_section<Scalar> bilinear_low_pass_section(Scalar omega, Scalar g1, Scalar g0)
{
    auto const beta1 = omega * g1;
    auto const beta0 = omega * omega * g0;
    auto const den0 = Scalar{1} + beta1 + beta0;
    auto const b0 = beta0 / den0;
    auto const pole_offset = Scalar{2} * beta1 / den0;
    return offset_section<Scalar>{
            .coeffs = biquad_coeffs<Scalar>{b0, Scalar{2} * b0, b0, Scalar{-2} + (Scalar{2} * beta1 + Scalar{4} * beta0) / den0, Scalar{1} - pole_offset},
            .dc_denominator = Scalar{4} * beta0 / den0,
            .pole_offset = pole_offset,
    };
}

// The RBJ notch rewritten exactly in Omega: the low-pass denominator with g1 = 1 / q and g0 = 1, b1 = a1, and
// b0 = b2 = (1 + Omega^2) / den0 stored as the offset 1 - beta1 / den0, so the numerator gap equals the design gap.
template <typename Scalar>
offset_section<Scalar> bilinear_notch_section(Scalar omega, Scalar q)
{
    auto section = bilinear_low_pass_section(omega, Scalar{1} / q, Scalar{1});
    auto const b0 = Scalar{1} - section.pole_offset / Scalar{2};
    section.coeffs.b0 = b0;
    section.coeffs.b1 = section.coeffs.a1;
    section.coeffs.b2 = b0;
    return section;
}

// The bilinear map of wc s / (s + wc) with k = 2 fs and wc = k Omega, divided through by k (1 + Omega): a1 is stored
// as the offset -1 + D1 from its limit, D1 = 1 + a1 = 2 Omega / (1 + Omega).
template <typename Scalar>
offset_section<Scalar> bilinear_differentiator_section(Scalar omega, Scalar k)
{
    auto const share = omega / (Scalar{1} + omega);
    auto const b0 = k * share;
    auto const gap = Scalar{2} * share;
    return offset_section<Scalar>{
            .coeffs = biquad_coeffs<Scalar>{b0, -b0, Scalar{0}, Scalar{-1} + gap, Scalar{0}},
            .dc_denominator = gap,
            .pole_offset = gap,
    };
}

// A section whose design gap is not finite cannot be judged against the floor, so finiteness is refused first.
template <typename Scalar>
ctrlpp::expected<biquad_coeffs<Scalar>, dsp_error> resolved_coefficients(offset_section<Scalar> const& section, Scalar resolution_factor)
{
    if(!coefficients_finite(section.coeffs) || !std::isfinite(section.dc_denominator))
        return ctrlpp::unexpected(dsp_error::non_finite_input);
    if(auto const resolved = check_coefficient_resolution(section.dc_denominator, resolution_factor); !resolved.has_value())
        return ctrlpp::unexpected(resolved.error());
    return section.coeffs;
}

// The unit-scale bilinear section of the prototype pole pair sigma +- j imag, g1 = -2 sigma and g0 = sigma^2 + imag^2,
// which is the bilinear map divided through by (2 fs)^2.
// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 7 (Chebyshev Type I analog prototype, bilinear pre-warp)
template <typename Scalar, std::size_t Order>
ctrlpp::expected<biquad_coeffs<Scalar>, dsp_error> chebyshev1_section(std::size_t k, Scalar sinh_v, Scalar cosh_v, Scalar omega, Scalar resolution_factor)
{
    auto const theta = std::numbers::pi_v<Scalar> * static_cast<Scalar>(2 * k + 1) / static_cast<Scalar>(2 * Order);
    auto const sigma = -sinh_v * std::sin(theta);
    auto const imag = cosh_v * std::cos(theta);
    return resolved_coefficients(bilinear_low_pass_section(omega, Scalar{-2} * sigma, sigma * sigma + imag * imag), resolution_factor);
}

}

#endif
