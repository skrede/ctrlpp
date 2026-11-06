#ifndef HPP_GUARD_CTRLPP_DSP_BIQUAD_H
#define HPP_GUARD_CTRLPP_DSP_BIQUAD_H

// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Sec. 6.4, Eq. (6.44) (transposed direct form II)
// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae for Audio EQ Biquad Filter Coefficients", 2005

#include "ctrlpp/dsp/dsp_types.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/dsp/detail/biquad_design.h"
#include "ctrlpp/dsp/detail/biquad_sections.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <numbers>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar>
class biquad
{
public:
    using scalar_type = Scalar;

    constexpr biquad()
        : m_c{}
        , m_w{}
    {
    }

    explicit constexpr biquad(biquad_coeffs<Scalar> c)
        : m_c{c}
        , m_w{}
    {
    }

    Scalar process(Scalar x)
    {
        auto const y = m_c.b0 * x + m_w[0];
        m_w[0] = m_c.b1 * x - m_c.a1 * y + m_w[1];
        m_w[1] = m_c.b2 * x - m_c.a2 * y;
        return y;
    }

    void reset() { m_w.fill(Scalar{0}); }

    void reset(Scalar value)
    {
        auto const denom = Scalar{1} + m_c.a1 + m_c.a2;
        auto const denom_scale = Scalar{1} + std::abs(m_c.a1) + std::abs(m_c.a2);
        if(std::abs(denom) < detail::biquad_singular_tol(denom_scale))
        {
            m_w.fill(Scalar{0});
            return;
        }
        auto const dc_gain = (m_c.b0 + m_c.b1 + m_c.b2) / denom;
        auto const y_ss = value * dc_gain;
        m_w[0] = value * (m_c.b1 + m_c.b2) - y_ss * (m_c.a1 + m_c.a2);
        m_w[1] = value * m_c.b2 - y_ss * m_c.a2;
    }

    biquad_coeffs<Scalar> const& coefficients() const { return m_c; }

    // @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (LPF section)
    // @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Sec. 7.3 (bilinear transform of analog prototypes, pre-warped by Eq. (7.26))
    static expected<biquad, dsp_error> low_pass(Scalar cutoff_hz, Scalar sample_hz, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
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

    // @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae", 2005 (notch section)
    // @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Sec. 7.2.2, Eqs. (7.18) and (7.26) (bilinear transform with pre-warping)
    static expected<biquad, dsp_error> notch(Scalar freq_hz, Scalar sample_hz, Scalar q, Scalar resolution_factor = detail::notch_resolution_count<Scalar>())
    {
        if(!std::isfinite(q))
            return unexpected(dsp_error::non_finite_input);
        if(auto const validated = detail::validate_biquad_design(freq_hz, sample_hz); !validated.has_value())
            return unexpected(validated.error());
        if(q <= Scalar{0})
            return unexpected(dsp_error::non_positive_q);
        if(auto const factor = detail::validate_resolution_factor(resolution_factor); !factor.has_value())
            return unexpected(factor.error());
        auto const omega = detail::prewarped_frequency(freq_hz, sample_hz);
        if(!omega.has_value())
            return unexpected(omega.error());
        auto const section = detail::resolved_coefficients(detail::bilinear_notch_section(*omega, q), resolution_factor);
        if(!section.has_value())
            return unexpected(section.error());
        return biquad{*section};
    }

    // @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Sec. 7.2.2, Eqs. (7.18) and (7.26) (bilinear transform with pre-warping)
    static expected<biquad, dsp_error> dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz, Scalar resolution_factor = detail::first_order_resolution_count<Scalar>())
    {
        auto const omega = detail::validated_prewarp(bandwidth_hz, sample_hz, resolution_factor);
        if(!omega.has_value())
            return unexpected(omega.error());
        auto const k = Scalar{2} * sample_hz;
        auto const section = detail::resolved_coefficients(detail::bilinear_differentiator_section(*omega, k), resolution_factor);
        if(!section.has_value())
            return unexpected(section.error());
        return biquad{*section};
    }

private:
    biquad_coeffs<Scalar> m_c;
    std::array<Scalar, 2> m_w;
};

template <typename Scalar>
biquad(biquad_coeffs<Scalar>) -> biquad<Scalar>;

namespace detail
{

static_assert(discrete_filter<biquad<double>>);
static_assert(discrete_filter<biquad<float>>);

}

}

// The cascade is built from biquad and so follows it; including biquad.h brings the cascade with it.
#include "ctrlpp/dsp/cascaded_biquad.h"

#endif
