// What the oracles in this file decide.
//
// A filter is a fixed arithmetic expression over its own coefficients, so where
// the coefficients make the answer exact, the answer is asserted exactly:
//
//  * All-zero coefficients, a unit direct coefficient with everything else zero,
//    and a single unit tap all produce their outputs BITWISE. Each of those case
//    names already said "exactly" or "identity", and a tolerance said less than
//    the arithmetic guarantees.
//  * Divergence is asserted against the CLOSED FORM of the recursion, the
//    particular solution plus the two homogeneous modes the poles fix, not
//    against a round magnitude nine decades below the realized value.
//  * Boundedness is asserted against the filter's own bounded-input gain -- the
//    absolute sum of the impulse response of the difference equation its
//    realized coefficients define, formed independently in the test -- times the
//    input bound. That is the supremum a bounded input can produce, so it is the
//    bound, and no round number is involved.
//  * A poisoned sample LATCHES in a recursive filter and FLUSHES from a
//    finite-impulse one after its delay line has run out. Both are asserted, and
//    each case names the other, because the distinction is a real property of
//    the two families and neither case said it alone.
//  * A design sweep over in-domain specifications asserts that each one designs.
//    A body inside a presence test asserts nothing on the absent branch, so if
//    every specification had failed the case would have passed silently.
//
// What they deliberately do not decide. The Chebyshev design case compares
// against reference coefficients with a counted rounding margin and its own
// comment already disclaims that the margin is a tolerance on the design; it is
// left as it stands and cited here as the file-local template for the counted
// budgets above.

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

using Catch::Matchers::WithinRel;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

/// Absolute sum of the impulse response of the direct-form difference equation
/// the coefficients define.
///
/// This is the exact supremum of the output magnitude over all inputs bounded by
/// one, so it is the bound a "stays bounded" case needs; a peak frequency
/// response would be a lower bound on it and a round number is neither.
///
/// The recursion is written here from the coefficients rather than taken from
/// the filter under test, so the bound does not inherit a defect in the state
/// update it is being used to check.
///
/// The response is summed over `samples` terms and the tail is discarded. That
/// is not an approximation at this scale: both poles of a stable biquad have
/// magnitude the square root of the trailing denominator coefficient, 0.64 for
/// the design below, and its two hundred and fifty sixth power is below the
/// smallest normal double, so nothing beyond the window can reach the sum.
auto bounded_input_gain(ctrlpp::biquad_coeffs<double> const& c) -> double
{
    constexpr int samples = 256;
    double y1 = 0.0;
    double y2 = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double gain = 0.0;

    for(int n = 0; n < samples; ++n)
    {
        double const x = (n == 0) ? 1.0 : 0.0;
        double const y = c.b0 * x + c.b1 * x1 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
        gain += std::abs(y);
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
    }
    return gain;
}

}

// ── Biquad hardening ───────────────────────────────────────────────────────────

TEST_CASE("Biquad with all-zero coefficients", "[biquad][hardening][negative]")
{
    ctrlpp::biquad_coeffs<double> c{};
    ctrlpp::biquad<double> filter(c);

    // Every coefficient is zero, so every product is an exact product of zero and
    // the output is exactly zero. Nothing here can round.
    double y = filter.process(1.0);
    REQUIRE(y == 0.0);

    y = filter.process(5.0);
    REQUIRE(y == 0.0);
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

    // Subsequent outputs should also be NaN (state contaminated). The poison
    // LATCHES here: the recursion feeds each output back into the state, so no
    // number of clean samples flushes it and only a reset recovers. The
    // finite-impulse case below is the counterpart, where the poison leaves with
    // the delay line instead, and the two together are what make the distinction
    // between the families a tested property rather than folklore.
    y = lp->process(1.0);
    CHECK(std::isnan(y));

    for(int k = 0; k < 8; ++k)
    {
        CAPTURE(k);
        CHECK(std::isnan(lp->process(1.0)));
    }
}

TEST_CASE("Biquad with unstable poles diverges", "[biquad][hardening][negative]")
{
    // Poles at z = +2.0 and z = +1.5 -- well outside unit circle
    // Denominator: z^2 - 3.5z + 3 => a1 = -3.5, a2 = 3.0
    ctrlpp::biquad_coeffs<double> c{.b0 = 1.0, .b1 = 0.0, .b2 = 0.0, .a1 = -3.5, .a2 = 3.0};
    ctrlpp::biquad<double> filter(c);

    // Feed constant input; output magnitude should grow
    constexpr int steps = 50;
    double y = 0.0;
    for (int i = 0; i < steps; ++i) {
        y = filter.process(1.0);
    }

    // The response has a closed form, so the case asserts it rather than a round
    // magnitude. The recursion is y_n = x_n + 3.5 y_{n-1} - 3 y_{n-2}; with a
    // constant unit input the particular solution is 1 / (1 - 3.5 + 3) = 2, and
    // the homogeneous part is a combination of the two pole powers. Matching the
    // first two outputs from rest, y_0 = 1 and y_1 = 4.5, fixes the combination:
    //     y_n = 2 + 8 * 2^n - 9 * 1.5^n.
    // At n = 49 that is 4.5036e15, so the bound this replaces sat nine decades
    // below the value it was checking and passed for a filter that barely moved.
    //
    // The budget is relative because the arithmetic and the answer grow together
    // in an unstable recursion: eight rounded operations per sample along the
    // three assignments of the transposed form, over fifty samples, plus six for
    // the closed form -- two power calls worth one unit in the last place each,
    // two products and two sums. Measured: half a unit in the last place.
    constexpr int recursion_ops_per_sample = 8;
    constexpr int closed_form_ops = 6;
    constexpr int divergence_ops = recursion_ops_per_sample * steps + closed_form_ops;

    double const expected = 2.0 + 8.0 * std::pow(2.0, steps - 1) - 9.0 * std::pow(1.5, steps - 1);
    CAPTURE(y, expected);
    CHECK_THAT(y, WithinRel(expected, divergence_ops * eps));
}

TEST_CASE("Unity gain biquad passes input through exactly", "[biquad][hardening][precision]")
{
    ctrlpp::biquad_coeffs<double> c{.b0 = 1.0, .b1 = 0.0, .b2 = 0.0, .a1 = 0.0, .a2 = 0.0};
    ctrlpp::biquad<double> filter(c);

    // The direct coefficient is one and every other is zero, so the output is the
    // input with no arithmetic that can round. The case name said "exactly" and
    // the assertion then allowed slack; it no longer does.
    REQUIRE(filter.process(3.14) == 3.14);
    REQUIRE(filter.process(-2.0) == -2.0);
    REQUIRE(filter.process(0.0) == 0.0);
}

TEST_CASE("Stable biquad output stays bounded over many samples", "[biquad][hardening][stability]")
{
    auto lp = ctrlpp::biquad<double>::low_pass(100.0, 1000.0);
    REQUIRE(lp.has_value());

    constexpr double input_bound = 1.0;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-input_bound, input_bound);

    double max_output = 0.0;
    for (int i = 0; i < 10000; ++i) {
        double y = lp->process(dist(gen));
        max_output = std::max(max_output, std::abs(y));
    }

    // The bound is the filter's own bounded-input gain, computed from the
    // coefficients it actually realized, times the bound on the input. That is
    // the exact supremum over every admissible input, so it is the tightest
    // statement "stays bounded" can make: 1.105 here, against the 10.0 it
    // replaces, with the realized peak at 0.826. A separate finiteness assertion
    // adds nothing once a finite bound is asserted, so it is gone.
    double const gain_bound = bounded_input_gain(lp->coefficients());
    CAPTURE(max_output, gain_bound);
    REQUIRE(max_output <= gain_bound * input_bound);
}

// ── FIR hardening ──────────────────────────────────────────────────────────────

TEST_CASE("FIR with all-zero coefficients", "[fir][hardening][negative]")
{
    ctrlpp::fir<double, 4> filter(std::array<double, 4>{0.0, 0.0, 0.0, 0.0});

    // Exact zeros, as above: every tap contributes an exact product of zero.
    REQUIRE(filter.process(1.0) == 0.0);
    REQUIRE(filter.process(5.0) == 0.0);
}

TEST_CASE("FIR poison flushes out of the delay line", "[fir][hardening][negative]")
{
    constexpr std::size_t taps = 3;
    ctrlpp::fir<double, taps> filter(std::array<double, taps>{0.25, 0.5, 0.25});

    filter.process(1.0);
    double y = filter.process(std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(y));

    // The half the case was missing, and the reason it matters: a finite-impulse
    // filter has no feedback, so the poisoned sample occupies exactly as many
    // outputs as it occupies taps and then leaves. Three taps, so the poisoned
    // output and the two after it are lost and the fourth is clean again. The
    // recursive case above is the counterpart, where the same injection never
    // leaves. Nothing tested that distinction before, and it is the operational
    // difference between the two families under a bad sample.
    for(std::size_t k = 1; k < taps; ++k)
    {
        CAPTURE(k);
        CHECK(std::isnan(filter.process(1.0)));
    }

    // Clean, and exactly so: the taps sum to one and each is a negative power of
    // two, so a unit input produces exactly one with no rounding anywhere.
    CHECK(filter.process(1.0) == 1.0);
    CHECK(filter.process(1.0) == 1.0);
}

TEST_CASE("Single-tap FIR with coefficient 1.0 is identity", "[fir][hardening][precision]")
{
    ctrlpp::fir<double, 1> filter(std::array<double, 1>{1.0});

    // One tap of one: the output is the input, bit for bit, and "identity" in the
    // case name is a claim about equality rather than about closeness.
    REQUIRE(filter.process(3.14) == 3.14);
    REQUIRE(filter.process(-7.0) == -7.0);
    REQUIRE(filter.process(0.0) == 0.0);
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
        // Every value in the list is strictly inside the design domain, so every
        // one of them must design. Skipping the failures asserted nothing on the
        // branch that skipped: had all six failed, the case would have reported
        // success having checked no coefficient at all.
        REQUIRE(design.has_value());

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

TEST_CASE("Biquad reset with a degenerate denominator zeroes the state",
          "[biquad][hardening][coverage]")
{
    // 1 + a1 + a2 is exactly zero here, so the steady state the reset would place
    // the filter in does not exist and the documented fallback is to zero the
    // state instead.
    ctrlpp::biquad_coeffs<double> c{.b0 = 1.0, .b1 = 0.0, .b2 = 0.0, .a1 = -0.5, .a2 = -0.5};
    ctrlpp::biquad<double> filter(c);

    filter.process(1.0);
    filter.reset(5.0);

    // The contract the case's own comment already stated and the assertion then
    // did not check. With a zeroed state and a unit direct coefficient, a zero
    // input produces exactly zero: the first output reads the first state word
    // and the second reads the second, so two consecutive exact zeros are what
    // proves BOTH words were cleared. Finiteness proved neither, and named no
    // property at all.
    REQUIRE(filter.process(0.0) == 0.0);
    REQUIRE(filter.process(0.0) == 0.0);
}
