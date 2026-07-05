// This anchor asserts the closed-form frequency response of two `biquad`
// factories against their textbook targets, rather than the loose amplitude
// thresholds used elsewhere in the suite.
//
// A digital biquad's magnitude response is evaluated exactly (up to floating
// point rounding) from its own coefficients by substituting z = e^{jOmega}
// into H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2), where
// Omega = 2*pi*f/fs is the digital frequency (Oppenheim &amp; Schafer,
// "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 5). This gives a
// purely mathematical oracle independent of the filter's own internal
// process/update path.
//
// A maximally-flat (Butterworth) low-pass section has a textbook property
// at its cutoff: |H(fc)| = 1/sqrt(2), i.e. -3.0103 dB, and the magnitude is
// monotonically non-increasing from DC up to the cutoff (no passband
// peaking) -- Oppenheim &amp; Schafer, Ch. 7, and Bristow-Johnson, "Cookbook
// Formulae for Audio EQ Biquad Filter Coefficients", 2005. The measured
// response here shows a passband peak above 0 dB and a +3 dB (not -3 dB)
// value at the nominal cutoff, so this section is held green with
// `[!shouldfail]` until the underlying quality-factor selection is
// corrected.
//
// A "dirty derivative" (a differentiator low-pass filtered above its own
// bandwidth) behaves, at frequencies well below that bandwidth, like an
// ideal differentiator: |H(f)| ~ 2*pi*f. This is the textbook low-frequency
// limit of a first-order high-pass-shaped differentiator, not a property
// derived from this codebase. The measured response here is off by roughly
// two orders of magnitude from that target, so this section is also held
// green with `[!shouldfail]` until corrected.

#include "ctrlpp/dsp/biquad.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <complex>
#include <numbers>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

/// Exact digital magnitude |H(e^{jOmega})| at Omega = 2*pi*f/fs, evaluated
/// directly from the biquad's own reported coefficients.
auto magnitude_linear(const biquad_coeffs<double>& c, double f, double fs) -> double
{
    const double omega = 2.0 * std::numbers::pi * f / fs;
    const std::complex<double> z_inv = std::exp(std::complex<double>(0.0, -omega));
    const std::complex<double> num = c.b0 + c.b1 * z_inv + c.b2 * z_inv * z_inv;
    const std::complex<double> den = 1.0 + c.a1 * z_inv + c.a2 * z_inv * z_inv;
    return std::abs(num / den);
}

auto magnitude_db(const biquad_coeffs<double>& c, double f, double fs) -> double
{
    return 20.0 * std::log10(magnitude_linear(c, f, fs));
}

} // namespace

TEST_CASE("Butterworth low-pass biquad reaches -3.01 dB at cutoff with a monotone passband", "[dsp][anchor][!shouldfail]")
{
    const double fc = 100.0;
    const double fs = 1000.0;
    const auto lp = biquad<double>::low_pass(fc, fs);
    const auto c = lp.coefficients();

    // Rounding-op margin: evaluating |H(e^{jOmega})| from coefficients chains
    // a handful of complex multiply-adds and one complex division; each
    // contributes up to one ULP of rounding at the magnitude's own O(1)
    // scale. Converting that linear-magnitude rounding into a decibel
    // tolerance scales it by the local derivative of 20*log10(x)/ln(10) at
    // the target magnitude 1/sqrt(2).
    constexpr double rounding_op_margin = 16.0;
    const double eps = std::numeric_limits<double>::epsilon();
    const double target_mag = 1.0 / std::numbers::sqrt2;
    const double db_per_unit_mag = 20.0 / (std::log(10.0) * target_mag);
    const double tol_db = db_per_unit_mag * rounding_op_margin * eps;

    // -10*log10(2) is the exact value of 20*log10(1/sqrt(2)): the textbook
    // half-power point, derived here rather than typed as a rounded literal.
    const double target_db = -10.0 * std::log10(2.0);
    const double measured_db = magnitude_db(c, fc, fs);
    CAPTURE(measured_db, target_db);
    REQUIRE_THAT(measured_db, WithinAbs(target_db, tol_db));

    // Passband monotonicity: from just above DC up to the cutoff, the
    // magnitude of a maximally-flat low-pass must never increase.
    constexpr int num_steps = 200;
    const double tol_lin = rounding_op_margin * eps;
    double prev_mag = magnitude_linear(c, fc / num_steps, fs);
    for(int i = 2; i <= num_steps; ++i)
    {
        const double f = fc * static_cast<double>(i) / static_cast<double>(num_steps);
        const double mag = magnitude_linear(c, f, fs);
        CAPTURE(f, mag, prev_mag);
        REQUIRE(mag <= prev_mag + tol_lin);
        prev_mag = mag;
    }
}

TEST_CASE("dirty_derivative magnitude tracks 2*pi*f well below its bandwidth", "[dsp][anchor][!shouldfail]")
{
    const double bandwidth_hz = 50.0;
    const double fs = 1000.0;
    const auto dd = biquad<double>::dirty_derivative(bandwidth_hz, fs);
    const auto c = dd.coefficients();

    // The bilinear transform maps the digital frequency Omega = 2*pi*f/fs to
    // an exact pre-warped analog frequency omega_a = 2*fs*tan(Omega/2)
    // (Oppenheim &amp; Schafer, Ch. 7): a bug-free digital biquad evaluates
    // H(e^{jOmega}) exactly equal to its continuous-time prototype evaluated
    // at omega_a. The dirty-derivative prototype is a differentiator rolled
    // off above a cutoff wc: H_analog(jw) = wc*jw/(jw+wc), so
    // |H_analog(jw)| = w/sqrt(1+(w/wc)^2). Two purely mathematical
    // inequalities -- neither derived from the biquad code under test --
    // bound how far the true digital magnitude may legitimately differ from
    // the textbook low-frequency differentiator limit 2*pi*f:
    //   (a) 1/sqrt(1+y) >= 1 - y/2 for y >= 0 (the function
    //       1/sqrt(1+y) - 1 + y/2 is zero at y=0 with non-negative
    //       derivative for y >= 0), applied at y = (omega_a/wc)^2, which
    //       bounds the differentiator's own rolloff below its cutoff;
    //   (b) tan(x) <= x + x^3 for x in [0,1] (the discarded Taylor terms
    //       x^3/3 + 2x^5/15 + ... are dominated by x^3 in this range),
    //       applied at x = pi*f/fs, which bounds the pre-warping frequency
    //       shift between the digital and analog frequency axes.
    const double wc = 2.0 * std::numbers::pi * bandwidth_hz;
    const double eps = std::numeric_limits<double>::epsilon();
    constexpr double rounding_op_margin = 16.0;

    for(double f : {1.0, 2.0, 3.0, 5.0})
    {
        const double x = std::numbers::pi * f / fs;
        const double omega_a = 2.0 * fs * std::tan(x);
        const double y = (omega_a / wc) * (omega_a / wc);
        const double rel_tol_rolloff = 0.5 * y;
        const double rel_tol_prewarp = x * x;
        const double rel_tol = rel_tol_rolloff + rel_tol_prewarp + rounding_op_margin * eps;

        const double expected = 2.0 * std::numbers::pi * f;
        const double measured = magnitude_linear(c, f, fs);

        CAPTURE(f, measured, expected, rel_tol);
        REQUIRE(measured <= expected * (1.0 + rel_tol));
        REQUIRE(measured >= expected * (1.0 - rel_tol));
    }
}
