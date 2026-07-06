#include "hardening_helpers.h"

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/fir.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <random>

using Catch::Matchers::WithinAbs;

// ── Biquad hardening ───────────────────────────────────────────────────────────

TEST_CASE("Biquad with all-zero coefficients", "[biquad][hardening][negative]")
{
    ctrlpp::biquad_coeffs<double> c{};
    ctrlpp::biquad<double> filter(c);

    double y = filter.process(1.0);
    REQUIRE_THAT(y, WithinAbs(0.0, 1e-15));

    y = filter.process(5.0);
    REQUIRE_THAT(y, WithinAbs(0.0, 1e-15));
}

TEST_CASE("Biquad with NaN input sample", "[biquad][hardening][negative]")
{
    auto lp = ctrlpp::biquad<double>::low_pass(100.0, 1000.0);
    REQUIRE(lp.has_value());

    // Process some valid samples first
    for (int i = 0; i < 10; ++i) {
        lp->process(1.0);
    }

    // Inject NaN
    double y = lp->process(std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(y));

    // Subsequent outputs should also be NaN (state contaminated)
    y = lp->process(1.0);
    CHECK(std::isnan(y));
}

TEST_CASE("Biquad with unstable poles diverges", "[biquad][hardening][negative]")
{
    // Poles at z = +2.0 and z = +1.5 -- well outside unit circle
    // Denominator: z^2 - 3.5z + 3 => a1 = -3.5, a2 = 3.0
    ctrlpp::biquad_coeffs<double> c{.b0 = 1.0, .b1 = 0.0, .b2 = 0.0, .a1 = -3.5, .a2 = 3.0};
    ctrlpp::biquad<double> filter(c);

    // Feed constant input; output magnitude should grow
    double y = 0.0;
    for (int i = 0; i < 50; ++i) {
        y = filter.process(1.0);
    }
    // After 50 steps with poles outside unit circle, magnitude should be large
    CHECK(std::abs(y) > 1e6);
}

TEST_CASE("Unity gain biquad passes input through exactly", "[biquad][hardening][precision]")
{
    ctrlpp::biquad_coeffs<double> c{.b0 = 1.0, .b1 = 0.0, .b2 = 0.0, .a1 = 0.0, .a2 = 0.0};
    ctrlpp::biquad<double> filter(c);

    REQUIRE_THAT(filter.process(3.14), WithinAbs(3.14, 1e-15));
    REQUIRE_THAT(filter.process(-2.0), WithinAbs(-2.0, 1e-15));
    REQUIRE_THAT(filter.process(0.0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("Stable biquad output stays bounded over many samples", "[biquad][hardening][stability]")
{
    auto lp = ctrlpp::biquad<double>::low_pass(100.0, 1000.0);
    REQUIRE(lp.has_value());

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    double max_output = 0.0;
    for (int i = 0; i < 10000; ++i) {
        double y = lp->process(dist(gen));
        max_output = std::max(max_output, std::abs(y));
    }

    // Low-pass with bounded input should produce bounded output
    REQUIRE(max_output < 10.0);
    REQUIRE(std::isfinite(max_output));
}

// ── FIR hardening ──────────────────────────────────────────────────────────────

TEST_CASE("FIR with all-zero coefficients", "[fir][hardening][negative]")
{
    ctrlpp::fir<double, 4> filter(std::array<double, 4>{0.0, 0.0, 0.0, 0.0});

    REQUIRE_THAT(filter.process(1.0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(filter.process(5.0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("FIR with NaN input", "[fir][hardening][negative]")
{
    ctrlpp::fir<double, 3> filter(std::array<double, 3>{0.25, 0.5, 0.25});

    filter.process(1.0);
    double y = filter.process(std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(y));
}

TEST_CASE("Single-tap FIR with coefficient 1.0 is identity", "[fir][hardening][precision]")
{
    ctrlpp::fir<double, 1> filter(std::array<double, 1>{1.0});

    REQUIRE_THAT(filter.process(3.14), WithinAbs(3.14, 1e-15));
    REQUIRE_THAT(filter.process(-7.0), WithinAbs(-7.0, 1e-15));
    REQUIRE_THAT(filter.process(0.0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("make_butterworth rejects invalid designs with the specific error", "[biquad][hardening][error]")
{
    auto zero_fs = ctrlpp::make_butterworth<4>(100.0, 0.0);
    REQUIRE(!zero_fs.has_value());
    CHECK(zero_fs.error() == ctrlpp::dsp_error::non_positive_sample_rate);

    // The boundary itself is rejected: the design interval is open, (0, fs/2).
    auto at_nyquist = ctrlpp::make_butterworth<4>(500.0, 1000.0);
    REQUIRE(!at_nyquist.has_value());
    CHECK(at_nyquist.error() == ctrlpp::dsp_error::cutoff_exceeds_nyquist);

    auto nan_cutoff = ctrlpp::make_butterworth<4>(std::numeric_limits<double>::quiet_NaN(), 1000.0);
    REQUIRE(!nan_cutoff.has_value());
    CHECK(nan_cutoff.error() == ctrlpp::dsp_error::non_finite_input);
}

TEST_CASE("make_chebyshev1 rejects invalid designs with the specific error", "[biquad][hardening][error]")
{
    auto negative_fs = ctrlpp::make_chebyshev1<4>(100.0, -1000.0, 1.0);
    REQUIRE(!negative_fs.has_value());
    CHECK(negative_fs.error() == ctrlpp::dsp_error::non_positive_sample_rate);

    auto above_nyquist = ctrlpp::make_chebyshev1<4>(600.0, 1000.0, 1.0);
    REQUIRE(!above_nyquist.has_value());
    CHECK(above_nyquist.error() == ctrlpp::dsp_error::cutoff_exceeds_nyquist);

    auto nan_ripple = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, std::numeric_limits<double>::quiet_NaN());
    REQUIRE(!nan_ripple.has_value());
    CHECK(nan_ripple.error() == ctrlpp::dsp_error::non_finite_input);
}

TEST_CASE("Biquad reset with near-zero denominator", "[biquad][hardening][coverage]")
{
    // a1 + a2 + 1 ~ 0 => denominator is near zero
    ctrlpp::biquad_coeffs<double> c{.b0 = 1.0, .b1 = 0.0, .b2 = 0.0, .a1 = -0.5, .a2 = -0.5};
    ctrlpp::biquad<double> filter(c);

    filter.process(1.0);
    filter.reset(5.0); // Should handle near-zero denominator gracefully

    // After reset with degenerate denominator, state should be zeroed
    double y = filter.process(0.0);
    CHECK(std::isfinite(y));
}
