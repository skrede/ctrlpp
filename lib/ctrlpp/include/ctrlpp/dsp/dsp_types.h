#ifndef HPP_GUARD_CTRLPP_DSP_DSP_TYPES_H
#define HPP_GUARD_CTRLPP_DSP_DSP_TYPES_H

// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 4 (sampling and the Nyquist criterion)

namespace ctrlpp
{

enum class dsp_error
{
    non_positive_sample_rate,
    cutoff_exceeds_nyquist,
    non_positive_q,
    non_positive_ripple,
    non_finite_input,
    invalid_resolution_factor,
    cutoff_below_coefficient_resolution,
};

template <typename Scalar>
struct biquad_coeffs
{
    constexpr biquad_coeffs()
        : b0{}
        , b1{}
        , b2{}
        , a1{}
        , a2{}
    {
    }

    constexpr biquad_coeffs(Scalar n0, Scalar n1, Scalar n2, Scalar d1, Scalar d2)
        : b0{n0}
        , b1{n1}
        , b2{n2}
        , a1{d1}
        , a2{d2}
    {
    }

    Scalar b0, b1, b2;
    Scalar a1, a2;
};

}

#endif
