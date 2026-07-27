#include "hardening_helpers.h"

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/fir.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstddef>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

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

TEST_CASE("make_chebyshev1 rejects a non-positive ripple specification",
          "[biquad][hardening][error]")
{
    // Exact domain condition, no tolerance involved: the ripple factor is
    // eps = sqrt(10^(ripple_db / 10) - 1), whose radicand is non-positive for
    // every ripple_db <= 0, and at exactly zero the following asinh(1 / eps)
    // takes an infinite argument. A non-positive ripple is therefore outside
    // the design's domain, not merely inaccurate. It is distinct from a
    // non-finite specification, which keeps its own enumerator.
    SECTION("a ripple of exactly zero")
    {
        auto zero_ripple = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, 0.0);
        REQUIRE(!zero_ripple.has_value());
        CHECK(zero_ripple.error() == ctrlpp::dsp_error::non_positive_ripple);
    }

    SECTION("a negative ripple")
    {
        auto negative_ripple = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, -1.0);
        REQUIRE(!negative_ripple.has_value());
        CHECK(negative_ripple.error() == ctrlpp::dsp_error::non_positive_ripple);
    }

    SECTION("a negative zero ripple")
    {
        auto negative_zero_ripple = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, -0.0);
        REQUIRE(!negative_zero_ripple.has_value());
        CHECK(negative_zero_ripple.error() == ctrlpp::dsp_error::non_positive_ripple);
    }
}

TEST_CASE("make_chebyshev1 leaves a conforming design unchanged",
          "[biquad][hardening][precision]")
{
    // Rounding-op margin for the Chebyshev Type I design chain. Each reference
    // coefficient below terminates a fixed closed-form chain: the ripple factor
    // (pow, subtract, sqrt), the prototype parameter (reciprocal, asinh, divide,
    // sinh, cosh), the pre-warped cutoff (tan, two products), the pole
    // coordinates (two products and a three-term sum of squares), the bilinear
    // denominator (two products and a three-term sum), its reciprocal, the
    // section coefficient product, and the DC-gain normalization across both
    // sections (two four-term quotients, their product, the target gain, and a
    // final scaling). Counting one rounding per arithmetic operation and one
    // unit in the last place per library transcendental gives 34 roundings along
    // the longest chain; the margin is that count times the scalar's machine
    // epsilon, applied relative to the magnitude of the coefficient compared.
    // It is not a tolerance on the design: the arithmetic here is unchanged, so
    // the same toolchain reproduces these values exactly. The margin exists only
    // so that another library's one-unit transcendental differences cannot fail
    // a bit-exact comparison.
    constexpr double design_roundings = 34.0;
    constexpr double design_margin = design_roundings * std::numeric_limits<double>::epsilon();

    auto const design = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, 1.0);
    REQUIRE(design.has_value());

    // Reference coefficients of the 4th-order, 1 dB passband ripple design at
    // a 100 Hz cutoff and a 1000 Hz sample rate.
    constexpr std::array<ctrlpp::biquad_coeffs<double>, 2> reference{
        ctrlpp::biquad_coeffs<double>{
            .b0 = 0.077686820472834692,
            .b1 = 0.15537364094566938,
            .b2 = 0.077686820472834692,
            .a1 = -1.49955449681044,
            .a2 = 0.84821868171669568,
        },
        ctrlpp::biquad_coeffs<double>{
            .b0 = 0.023627564635016508,
            .b1 = 0.047255129270033017,
            .b2 = 0.023627564635016508,
            .a1 = -1.5547851795965146,
            .a2 = 0.64929543813658086,
        },
    };

    for(std::size_t k = 0; k < reference.size(); ++k)
    {
        CAPTURE(k);
        auto const& c = design->section(k).coefficients();
        auto const& r = reference[k];
        CHECK_THAT(c.b0, WithinRel(r.b0, design_margin));
        CHECK_THAT(c.b1, WithinRel(r.b1, design_margin));
        CHECK_THAT(c.b2, WithinRel(r.b2, design_margin));
        CHECK_THAT(c.a1, WithinRel(r.a1, design_margin));
        CHECK_THAT(c.a2, WithinRel(r.a2, design_margin));
    }
}

TEST_CASE("make_chebyshev1 rejects a design whose coefficients degenerate",
          "[biquad][hardening][error]")
{
    // The ripple guard is a domain bound on the specification; the coefficient
    // sweep before the success return is what catches the remaining parameter
    // combinations that reach a non-finite coefficient through a path the
    // domain bounds do not describe. Both cases below satisfy every input
    // domain check and still degenerate, so they reach the sweep, not a guard.
    SECTION("a positive ripple so small the ripple factor underflows to zero")
    {
        // 10^(ripple_db / 10) rounds to exactly one here, so the radicand
        // 10^(ripple_db / 10) - 1 is exactly zero, eps is zero, and the
        // following asinh(1 / eps) takes an infinite argument. The
        // specification is strictly positive, so the ripple guard passes it.
        auto const design = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, 1e-300);
        REQUIRE(!design.has_value());
        CHECK(design.error() == ctrlpp::dsp_error::non_finite_input);
    }

    SECTION("a sample rate large enough to overflow the bilinear pre-warp")
    {
        // The pre-warped cutoff is wc = 2 fs tan(pi cutoff / fs); at this scale
        // it overflows, the section numerator becomes infinite, the reciprocal
        // of the denominator becomes zero, and their product is NaN. The cutoff
        // still lies strictly inside (0, fs / 2), so the Nyquist guard passes it.
        auto const design = ctrlpp::make_chebyshev1<4>(4.9e299, 1e300, 1.0);
        REQUIRE(!design.has_value());
        CHECK(design.error() == ctrlpp::dsp_error::non_finite_input);
    }
}

TEST_CASE("A successful chebyshev1 design never carries a non-finite coefficient",
          "[biquad][hardening][error]")
{
    // Ripple specifications spanning three decades, from a ripple far below the
    // resolution of any realizable passband to one that places the prototype
    // poles close to the imaginary axis. Every design that reports success must
    // hand back coefficients a filter can actually run.
    constexpr std::array<double, 6> ripples_db{1e-3, 1e-2, 0.1, 1.0, 3.0, 20.0};

    for(double ripple_db : ripples_db)
    {
        CAPTURE(ripple_db);
        auto const design = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, ripple_db);
        if(!design.has_value())
            continue;

        for(std::size_t k = 0; k < 2; ++k)
        {
            CAPTURE(k);
            auto const& c = design->section(k).coefficients();
            CHECK(std::isfinite(c.b0));
            CHECK(std::isfinite(c.b1));
            CHECK(std::isfinite(c.b2));
            CHECK(std::isfinite(c.a1));
            CHECK(std::isfinite(c.a2));
        }
    }
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
