#ifndef HPP_GUARD_CTRLPP_DSP_DETAIL_BIQUAD_DESIGN_H
#define HPP_GUARD_CTRLPP_DSP_DETAIL_BIQUAD_DESIGN_H

#include "ctrlpp/dsp/dsp_types.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <numbers>

namespace ctrlpp::detail
{

// Scaled by Scalar's own epsilon and the magnitude of the summed operands, so it holds at float and double alike. The
// default coefficient is twice the number of terms in the tested sum 1 + a1 + a2.
template <typename Scalar>
Scalar biquad_singular_tol(Scalar scale, Scalar coeff = Scalar{6})
{
    return coeff * std::numeric_limits<Scalar>::epsilon() * scale;
}

template <typename Scalar>
ctrlpp::expected<void, dsp_error> validate_biquad_design(Scalar freq_hz, Scalar sample_hz)
{
    if(!std::isfinite(freq_hz) || !std::isfinite(sample_hz))
        return ctrlpp::unexpected(dsp_error::non_finite_input);
    if(sample_hz <= Scalar{0})
        return ctrlpp::unexpected(dsp_error::non_positive_sample_rate);
    if(freq_hz <= Scalar{0} || freq_hz >= sample_hz / Scalar{2})
        return ctrlpp::unexpected(dsp_error::cutoff_exceeds_nyquist);
    return {};
}

// The representation error of a section's 1 + a1 + a2 with a1 stored near -2 and a2 near 1: half an ulp of a1 in
// [1, 2) is eps / 2 and half an ulp of a2 in [1/2, 1) is eps / 4, 3/4 eps in all (oppenheim2010dsp Sec. 6.8,
// coefficient quantization with poles near z = 1).
template <typename Scalar>
Scalar second_order_resolution_count()
{
    return Scalar{3} / Scalar{4};
}

// The representation error of a notch section's DC gain, the ratio of its numerator gap b0 + b1 + b2 to its design gap
// 1 + a1 + a2: a1 (eps / 2) and a2 (eps / 4) as in the low-pass section, b1 stored like a1 (eps / 2), and b0 and b2
// stored in [1/2, 1) (eps / 4 each), 7/4 eps in all.
template <typename Scalar>
Scalar notch_resolution_count()
{
    return Scalar{7} / Scalar{4};
}

// The representation error of a first-order section's design gap 1 + a1, with a1 stored in [-1, -1/2): half an ulp,
// eps / 4, the gap's only cancellation.
template <typename Scalar>
Scalar first_order_resolution_count()
{
    return Scalar{1} / Scalar{4};
}

template <typename Scalar>
ctrlpp::expected<void, dsp_error> validate_resolution_factor(Scalar resolution_factor)
{
    if(!std::isfinite(resolution_factor) || !(resolution_factor > Scalar{0}))
        return ctrlpp::unexpected(dsp_error::invalid_resolution_factor);
    return {};
}

template <typename Scalar>
ctrlpp::expected<void, dsp_error> check_coefficient_resolution(Scalar dc_denominator, Scalar resolution_factor)
{
    if(!(resolution_factor * std::numeric_limits<Scalar>::epsilon() < dc_denominator))
        return ctrlpp::unexpected(dsp_error::cutoff_below_coefficient_resolution);
    return {};
}

// Omega = tan(pi fc / fs), with the ratio formed first so pi times a large frequency cannot overflow. A result that
// is not finite and strictly positive means the frequency rounded onto zero or onto Nyquist in Scalar arithmetic.
template <typename Scalar>
ctrlpp::expected<Scalar, dsp_error> prewarped_frequency(Scalar freq_hz, Scalar sample_hz)
{
    auto const omega = std::tan(std::numbers::pi_v<Scalar> * (freq_hz / sample_hz));
    if(!std::isfinite(omega) || !(omega > Scalar{0}))
        return ctrlpp::unexpected(dsp_error::cutoff_exceeds_nyquist);
    return omega;
}

template <typename Scalar>
bool coefficients_finite(biquad_coeffs<Scalar> const& c)
{
    return std::isfinite(c.b0) && std::isfinite(c.b1) && std::isfinite(c.b2) && std::isfinite(c.a1) && std::isfinite(c.a2);
}

template <typename Scalar>
ctrlpp::expected<Scalar, dsp_error> validated_prewarp(Scalar freq_hz, Scalar sample_hz, Scalar resolution_factor)
{
    if(auto const validated = validate_biquad_design(freq_hz, sample_hz); !validated.has_value())
        return ctrlpp::unexpected(validated.error());
    if(auto const factor = validate_resolution_factor(resolution_factor); !factor.has_value())
        return ctrlpp::unexpected(factor.error());
    return prewarped_frequency(freq_hz, sample_hz);
}

}

#endif
