#ifndef HPP_GUARD_CTRLPP_DSP_DSP_TYPES_H
#define HPP_GUARD_CTRLPP_DSP_DSP_TYPES_H

/// @brief Public types for the dsp filter design factories.
///
/// `dsp_error` enumerates the structured failure modes the biquad-family
/// design factories can produce. Together with the designed filter they form
/// the `ctrlpp::expected<Filter, dsp_error>` contract of `biquad::low_pass`,
/// `biquad::notch`, `biquad::dirty_derivative`, `make_butterworth`,
/// `make_chebyshev1`, and their `vector_biquad` counterparts.
///
/// @cite oppenheim2010dsp -- Oppenheim &amp; Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 4 (sampling and the Nyquist criterion)

namespace ctrlpp
{

/// @brief Structured failure modes for the biquad-family design factories.
///
///  * non_positive_sample_rate : sample_hz <= 0; a sampled-data design requires
///                               a strictly positive sample rate.
///  * cutoff_exceeds_nyquist   : the design frequency (cutoff, notch center, or
///                               differentiator bandwidth) lies outside the open
///                               interval (0, sample_hz / 2). A discrete-time
///                               filter can only realize a response strictly
///                               below half the sample rate; at or above it the
///                               design frequency aliases (Nyquist criterion,
///                               oppenheim2010dsp Ch. 4).
///  * non_positive_q           : notch quality factor q <= 0; the design
///                               bandwidth alpha = sin(w0) / (2 q) requires
///                               q > 0.
///  * non_positive_ripple      : Chebyshev Type I passband ripple ripple_db <= 0.
///                               The ripple factor is
///                               eps = sqrt(10^(ripple_db / 10) - 1), whose
///                               radicand is non-positive for every
///                               ripple_db <= 0, and at exactly zero the
///                               following asinh(1 / eps) takes an infinite
///                               argument. An equiripple passband is defined by
///                               a strictly positive ripple, so this is an exact
///                               domain bound and carries no tolerance.
///  * non_finite_input         : a design parameter (frequency, sample rate,
///                               quality factor, or ripple) is NaN or infinite,
///                               or the design chain produced a non-finite
///                               coefficient.
enum class dsp_error
{
    non_positive_sample_rate,
    cutoff_exceeds_nyquist,
    non_positive_q,
    non_positive_ripple,
    non_finite_input,
};

}

#endif
