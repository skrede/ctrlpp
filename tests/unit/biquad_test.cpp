#include "ctrlpp/dsp/biquad.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <numbers>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace {

auto sine_amplitude_after_settling(auto& filter, double freq_hz, double sample_hz,
                                   int total_samples, int settle_samples) -> double
{
    double max_val = 0.0;
    for (int i = 0; i < total_samples; ++i) {
        double t = static_cast<double>(i) / sample_hz;
        double x = std::sin(2.0 * std::numbers::pi * freq_hz * t);
        double y = filter.process(x);
        if (i >= settle_samples) {
            max_val = std::max(max_val, std::abs(y));
        }
    }
    return max_val;
}

} // namespace

TEST_CASE("biquad low-pass passes DC signal unchanged", "[biquad]")
{
    auto lp = biquad<double>::low_pass(100.0, 1000.0);
    double y = 0.0;
    for (int i = 0; i < 200; ++i) {
        y = lp.process(1.0);
    }
    CHECK_THAT(y, WithinAbs(1.0, 0.01));
}

TEST_CASE("biquad low-pass attenuates above cutoff", "[biquad]")
{
    auto lp_pass = biquad<double>::low_pass(100.0, 1000.0);
    auto lp_stop = biquad<double>::low_pass(100.0, 1000.0);

    double passband = sine_amplitude_after_settling(lp_pass, 50.0, 1000.0, 2000, 500);
    double stopband = sine_amplitude_after_settling(lp_stop, 400.0, 1000.0, 2000, 500);

    CHECK(passband > 0.5);
    CHECK(stopband < 0.2);
}

TEST_CASE("biquad notch rejects target frequency", "[biquad]")
{
    auto n_reject = biquad<double>::notch(60.0, 1000.0, 10.0);
    auto n_pass = biquad<double>::notch(60.0, 1000.0, 10.0);

    double reject = sine_amplitude_after_settling(n_reject, 60.0, 1000.0, 2000, 500);
    double pass = sine_amplitude_after_settling(n_pass, 200.0, 1000.0, 2000, 500);

    CHECK(reject < 0.1);
    CHECK(pass > 0.5);
}

TEST_CASE("biquad dirty_derivative rejects DC and passes high frequency", "[biquad]")
{
    double fs = 1000.0;
    auto dd = biquad<double>::dirty_derivative(50.0, fs);

    // DC rejection: constant input should produce zero output after settling
    double y_dc = 0.0;
    for (int i = 0; i < 500; ++i) {
        y_dc = dd.process(5.0);
    }
    CHECK(std::abs(y_dc) < 0.01);

    // Step response: derivative of step is impulse, first output should be large
    dd.reset();
    double y_step = dd.process(1.0);
    CHECK(std::abs(y_step) > 0.5); // initial step response is significant

    // After many samples of constant input, output decays toward zero
    for (int i = 0; i < 500; ++i) {
        y_step = dd.process(1.0);
    }
    CHECK(std::abs(y_step) < 0.01);
}

TEST_CASE("biquad reset clears internal state", "[biquad]")
{
    auto lp = biquad<double>::low_pass(100.0, 1000.0);

    // Process some samples to build up state
    for (int i = 0; i < 100; ++i) {
        lp.process(1.0);
    }

    lp.reset();
    double y = lp.process(1.0);
    // After reset, output should equal b0 * input (no state carry)
    auto c = lp.coefficients();
    CHECK_THAT(y, WithinAbs(c.b0 * 1.0, 1e-12));
}

TEST_CASE("biquad reset with value initializes to steady state", "[biquad]")
{
    auto lp = biquad<double>::low_pass(100.0, 1000.0);
    lp.reset(5.0);
    double y = lp.process(5.0);

    // For a low-pass filter, DC gain is 1.0, so steady state output for 5.0 should be near 5.0
    CHECK_THAT(y, WithinAbs(5.0, 0.1));
}

TEST_CASE("cascaded_biquad butterworth order 4 attenuates more than single biquad", "[biquad]")
{
    auto single = biquad<double>::low_pass(100.0, 1000.0);
    auto cascade = make_butterworth<4>(100.0, 1000.0);

    double single_stop = sine_amplitude_after_settling(single, 400.0, 1000.0, 2000, 500);
    double cascade_stop = sine_amplitude_after_settling(cascade, 400.0, 1000.0, 2000, 500);

    CHECK(cascade_stop < single_stop);
}

TEST_CASE("make_chebyshev1 has steeper rolloff than butterworth", "[biquad]")
{
    auto butter = make_butterworth<4>(100.0, 1000.0);
    auto cheby = make_chebyshev1<4>(100.0, 1000.0, 1.0);

    // At a moderately-above-cutoff frequency, Chebyshev should attenuate more
    double butter_stop = sine_amplitude_after_settling(butter, 250.0, 1000.0, 3000, 1000);
    double cheby_stop = sine_amplitude_after_settling(cheby, 250.0, 1000.0, 3000, 1000);

    CHECK(cheby_stop < butter_stop);
}
