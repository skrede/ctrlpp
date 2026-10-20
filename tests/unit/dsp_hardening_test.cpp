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
//  * The Chebyshev precision cases compare against the EXACT design, computed at
//    60 significant digits by a standard-library generator, with one-sided
//    counted rounding margins. The first section's numerator, which carries the
//    DC normalization, is judged through the stored cascade's DC gain, evaluated
//    in double-double against the exact DC target, so the margin covers the
//    normalization's own rounding and not the offsets' formation error divided by
//    each gap.

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
#include <numbers>

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

// Double-double arithmetic, so a cascade's DC gain can be evaluated from its stored coefficients with a relative error
// near 2^-104: TwoSum, Fast2Sum and the fma-based TwoProduct (Dekker, "A floating-point technique for extending the
// available precision", Numer. Math. 18, 1971; Muller et al., "Handbook of Floating-Point Arithmetic", 2nd ed., 2018,
// Ch. 4 and 14), with the quotient refined by one Newton correction.
struct double_double
{
    double hi;
    double lo;
};

auto two_sum(double a, double b) -> double_double
{
    double const s = a + b;
    double const b_virtual = s - a;
    return {s, (a - (s - b_virtual)) + (b - b_virtual)};
}

auto fast_two_sum(double a, double b) -> double_double
{
    double const s = a + b;
    return {s, b - (s - a)};
}

auto two_product(double a, double b) -> double_double
{
    double const p = a * b;
    return {p, std::fma(a, b, -p)};
}

auto add(double_double x, double_double y) -> double_double
{
    auto const s = two_sum(x.hi, y.hi);
    auto const t = two_sum(x.lo, y.lo);
    auto const v = fast_two_sum(s.hi, s.lo + t.hi);
    return fast_two_sum(v.hi, t.lo + v.lo);
}

auto multiply(double_double x, double_double y) -> double_double
{
    auto const p = two_product(x.hi, y.hi);
    return fast_two_sum(p.hi, p.lo + (x.hi * y.lo + x.lo * y.hi));
}

auto divide(double_double x, double_double y) -> double_double
{
    double const q = x.hi / y.hi;
    auto const residual = add(x, multiply(y, double_double{-q, 0.0}));
    return fast_two_sum(q, residual.hi / y.hi);
}

template <typename Scalar, std::size_t N>
auto stored_dc_gain(ctrlpp::cascaded_biquad<Scalar, N> const& cascade) -> double_double
{
    double_double gain{1.0, 0.0};
    for(std::size_t k = 0; k < N; ++k)
    {
        auto const& c = cascade.section(k).coefficients();
        auto const numerator = add(add({static_cast<double>(c.b0), 0.0}, {static_cast<double>(c.b1), 0.0}), {static_cast<double>(c.b2), 0.0});
        auto const denominator = add(add({1.0, 0.0}, {static_cast<double>(c.a1), 0.0}), {static_cast<double>(c.a2), 0.0});
        gain = multiply(gain, divide(numerator, denominator));
    }
    return gain;
}

auto relative_deviation(double_double value, double_double exact) -> double
{
    return std::abs(add(value, {-exact.hi, -exact.lo}).hi / exact.hi);
}

struct exact_section
{
    double_double a1;
    double_double a2;
    double_double b0;
};

struct counted_margins
{
    double offset;
    double numerator;
    double dc_gain;
};

/// The first section's numerator carries the DC normalization, so it is judged through the stored cascade's DC gain,
/// which does not see the offsets' formation error: the normalization divided by the same stored gaps.
template <typename Scalar>
void check_against_exact_design(ctrlpp::cascaded_biquad<Scalar, 2> const& design, std::array<exact_section, 2> const& exact, double_double dc_target, counted_margins const& counts)
{
    auto const unit = static_cast<double>(std::numeric_limits<Scalar>::epsilon());
    auto const deviation = [unit](Scalar value, double_double reference) {
        return relative_deviation({static_cast<double>(value), 0.0}, reference) / unit;
    };
    for(std::size_t k = 0; k < exact.size(); ++k)
    {
        CAPTURE(k);
        auto const& c = design.section(k).coefficients();
        CHECK(c.b1 == Scalar{2} * c.b0);
        CHECK(c.b2 == c.b0);
        CHECK(deviation(c.a1, exact[k].a1) <= counts.offset);
        CHECK(deviation(c.a2, exact[k].a2) <= counts.offset);
        if(k > 0)
            CHECK(deviation(c.b0, exact[k].b0) <= counts.numerator);
    }
    CHECK(relative_deviation(stored_dc_gain(design), dc_target) / unit <= counts.dc_gain);
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

TEST_CASE("biquad low_pass refuses a single-precision cutoff below the coefficient resolution", "[biquad][hardening][error]")
{
    // At 0.03 Hz / 1000 Hz the design gap 1 + a1 + a2 is about 0.3 float eps, below the 3/4 eps the stored
    // a1 and a2 resolve, so no float filter carries this design.
    auto const design = ctrlpp::biquad<float>::low_pass(0.03f, 1000.f);
    REQUIRE(!design.has_value());
    CHECK(design.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);
}

TEST_CASE("biquad low_pass refuses a resolution factor that is not finite and positive", "[biquad][hardening][error]")
{
    for(double factor : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), 0.0, -0.0, -1.0})
    {
        CAPTURE(factor);
        auto const design = ctrlpp::biquad<double>::low_pass(100.0, 1000.0, factor);
        REQUIRE(!design.has_value());
        CHECK(design.error() == ctrlpp::dsp_error::invalid_resolution_factor);
    }
}

TEST_CASE("biquad low_pass resolution factor sets the refusal floor", "[biquad][hardening][error]")
{
    // At 1 Hz / 1000 Hz the float design gap is about 330 eps: accepted at the default 3/4, refused when the caller
    // asks the gap to be resolved to a relative accuracy of 1e-3, a factor of 3/4 / 1e-3 = 750.
    REQUIRE(ctrlpp::biquad<float>::low_pass(1.0f, 1000.f).has_value());
    auto const strict = ctrlpp::biquad<float>::low_pass(1.0f, 1000.f, 750.0f);
    REQUIRE(!strict.has_value());
    CHECK(strict.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);
}

TEST_CASE("biquad low_pass refuses a cutoff ratio that rounds to zero in the pre-warp", "[biquad][hardening][error]")
{
    // 1e-38 / 1e30 lies inside (0, fs / 2) but the ratio underflows to zero in float, so the pre-warped frequency is
    // zero: the design frequency rounded onto the lower edge of the open interval.
    auto const design = ctrlpp::biquad<float>::low_pass(1e-38f, 1e30f);
    REQUIRE(!design.has_value());
    CHECK(design.error() == ctrlpp::dsp_error::cutoff_exceeds_nyquist);
}

TEST_CASE("biquad notch refuses a single-precision center below the coefficient resolution", "[biquad][hardening][error]")
{
    // At 0.03 Hz / 1000 Hz and q = 0.7 the design gap is about 0.3 float eps, below the 7/4 eps the stored numerator
    // and denominator resolve between them.
    auto const design = ctrlpp::biquad<float>::notch(0.03f, 1000.f, 0.7f);
    REQUIRE(!design.has_value());
    CHECK(design.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);

    // At 1 Hz / 1000 Hz the gap is about 330 float eps: accepted at the default, refused at a relative accuracy of
    // 1e-3, a factor of 7/4 / 1e-3 = 1750.
    REQUIRE(ctrlpp::biquad<float>::notch(1.0f, 1000.f, 0.7f).has_value());
    auto const strict = ctrlpp::biquad<float>::notch(1.0f, 1000.f, 0.7f, 1750.0f);
    REQUIRE(!strict.has_value());
    CHECK(strict.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);
}

TEST_CASE("biquad notch and dirty_derivative refuse a resolution factor that is not finite and positive", "[biquad][hardening][error]")
{
    for(double factor : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), 0.0, -0.0, -1.0})
    {
        CAPTURE(factor);
        auto const notch = ctrlpp::biquad<double>::notch(60.0, 1000.0, 10.0, factor);
        REQUIRE(!notch.has_value());
        CHECK(notch.error() == ctrlpp::dsp_error::invalid_resolution_factor);
        auto const derivative = ctrlpp::biquad<double>::dirty_derivative(50.0, 1000.0, factor);
        REQUIRE(!derivative.has_value());
        CHECK(derivative.error() == ctrlpp::dsp_error::invalid_resolution_factor);
    }
}

TEST_CASE("biquad notch refuses a quality factor whose reciprocal overflows", "[biquad][hardening][error]")
{
    // The smallest subnormal q is positive, but 1 / q overflows and the design would carry NaN coefficients.
    auto const single = ctrlpp::biquad<float>::notch(100.f, 1000.f, std::numeric_limits<float>::denorm_min());
    REQUIRE(!single.has_value());
    CHECK(single.error() == ctrlpp::dsp_error::non_finite_input);
    auto const twice = ctrlpp::biquad<double>::notch(100.0, 1000.0, std::numeric_limits<double>::denorm_min());
    REQUIRE(!twice.has_value());
    CHECK(twice.error() == ctrlpp::dsp_error::non_finite_input);
}

TEST_CASE("biquad dirty_derivative refuses a bandwidth below its first-order floor", "[biquad][hardening][error]")
{
    // At 1e-6 Hz / 1000 Hz the gap 1 + a1 is about 2 pi 1e-9, below the eps / 4 a1 resolves near -1.
    auto const design = ctrlpp::biquad<float>::dirty_derivative(1e-6f, 1000.f);
    REQUIRE(!design.has_value());
    CHECK(design.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);
}

TEST_CASE("biquad dirty_derivative refuses a gain that overflows", "[biquad][hardening][error]")
{
    // The derivative gain 2 sample_hz is not representable, so b0 is infinite.
    auto const twice = ctrlpp::biquad<double>::dirty_derivative(100.0, 1e308);
    REQUIRE(!twice.has_value());
    CHECK(twice.error() == ctrlpp::dsp_error::non_finite_input);
    auto const single = ctrlpp::biquad<float>::dirty_derivative(100.f, 3e38f);
    REQUIRE(!single.has_value());
    CHECK(single.error() == ctrlpp::dsp_error::non_finite_input);
}

TEST_CASE("biquad notch and dirty_derivative refuse a ratio that rounds to zero in the pre-warp", "[biquad][hardening][error]")
{
    auto const notch = ctrlpp::biquad<float>::notch(1e-38f, 1e30f, 0.7f);
    REQUIRE(!notch.has_value());
    CHECK(notch.error() == ctrlpp::dsp_error::cutoff_exceeds_nyquist);
    auto const derivative = ctrlpp::biquad<float>::dirty_derivative(1e-38f, 1e30f);
    REQUIRE(!derivative.has_value());
    CHECK(derivative.error() == ctrlpp::dsp_error::cutoff_exceeds_nyquist);
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

TEST_CASE("Single-precision cascades below the coefficient-resolution floor are refused", "[biquad][hardening][error]")
{
    // Butterworth sections all have g0 = 1, so at 3e-5 of the sample rate each gap is about 0.3 float eps.
    auto const butterworth = ctrlpp::make_butterworth<4>(0.03f, 1000.f);
    REQUIRE(!butterworth.has_value());
    CHECK(butterworth.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);

    // At 7e-5 the 1 dB Chebyshev cascade's first section (pole radius squared 0.987) clears the floor with a gap of
    // 1.6 float eps, and its second (0.279) does not, at 0.45 eps: every section is checked, the smallest gap deciding.
    auto const chebyshev = ctrlpp::make_chebyshev1<4>(0.07f, 1000.f, 1.0f);
    REQUIRE(!chebyshev.has_value());
    CHECK(chebyshev.error() == ctrlpp::dsp_error::cutoff_below_coefficient_resolution);

    for(float factor : {std::numeric_limits<float>::quiet_NaN(), 0.0f})
    {
        CAPTURE(factor);
        auto const invalid_butterworth = ctrlpp::make_butterworth<4>(100.f, 1000.f, factor);
        REQUIRE(!invalid_butterworth.has_value());
        CHECK(invalid_butterworth.error() == ctrlpp::dsp_error::invalid_resolution_factor);
        auto const invalid_chebyshev = ctrlpp::make_chebyshev1<4>(100.f, 1000.f, 1.0f, factor);
        REQUIRE(!invalid_chebyshev.has_value());
        CHECK(invalid_chebyshev.error() == ctrlpp::dsp_error::invalid_resolution_factor);
    }
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

TEST_CASE("make_chebyshev1 matches the exact design within its counted rounding",
          "[biquad][hardening][precision]")
{
    // The reference is the exact 4th-order, 1 dB design at a 100 Hz cutoff and a 1000 Hz sample rate, printed by
    // validation/oracles/chebyshev1/oracle_chebyshev1.py as the double nearest each value and the double nearest the
    // remainder. The margins are one-sided counts of the formation's running error bound (Higham, "Accuracy and
    // Stability of Numerical Algorithms", 2nd ed., 2002, Sec. 3.3): one eps per operation and per library call, each
    // call carrying its argument's error through its slope |x f'(x) / f(x)|, power-of-two scalings and Sterbenz-exact
    // differences free, the first-order sum rounded up so the ceiling absorbs the second-order terms. The largest
    // offset count is 10.97 (a2 of the second section), so 11; the second section's unnormalized numerator 35.07, so
    // 36; the DC gain 10.5 for the normalization (per section 1.75 for b0 + b1 + b2, 1 for the gap since 1 + a1 is
    // Sterbenz-exact, 1 for the quotient; then the product, the correction and its product) plus 3.14 for the target
    // 1 / sqrt(1 + eps^2), so 14. The double-double evaluation adds near 2^-104.
    constexpr std::array<exact_section, 2> exact{
        exact_section{
            .a1 = {-0x1.7fe2cdb25ede1p+0, -0x1.5c458107b7f71p-56},
            .a2 = {0x1.b249b813a8ac9p-1, -0x1.f3cfa944b02f0p-56},
            .b0 = {0x1.6508395dd5e0dp-4, 0x1.5fd5ab672fb3ep-59},
        },
        exact_section{
            .a1 = {-0x1.8e0666caac2a1p+0, -0x1.78d46386e4c90p-54},
            .a2 = {0x1.4c7073a07a2a8p-1, 0x1.228a08a43bd25p-57},
            .b0 = {0x1.831d30590eb2bp-6, -0x1.48322725d4ebdp-60},
        },
    };
    constexpr double_double dc_target{0x1.c8520affa0a4bp-1, -0x1.2e5103f561a6fp-55};

    auto const design = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, 1.0);
    REQUIRE(design.has_value());
    check_against_exact_design(*design, exact, dc_target, counted_margins{.offset = 11.0, .numerator = 36.0, .dc_gain = 14.0});
}

TEST_CASE("make_chebyshev1 matches the exact small-ripple design within its counted rounding in single precision",
          "[biquad][hardening][precision]")
{
    // Formed as 10^(r / 10) - 1, the radicand at 1e-4 dB cancels to about three significant digits in float and the
    // coefficients land thousands of float epsilons from the exact design. Through expm1 they stay within their
    // counted bounds, counted as in the case above in float eps against the exact design at the float ripple
    // 0x1.a36e2ep-14 from the same generator: the offsets 306.55 of their own magnitude, so 307, since a2 of the
    // second section is 0.076 and its offset 0.924 weighs twelve times; the unnormalized numerator 40.62, so 41; the
    // DC gain 11.94 for the normalization (1 + a1 is not Sterbenz-exact here, |a1| < 1/2) plus 2.50 for the target,
    // so 15.
    constexpr std::array<exact_section, 2> exact{
        exact_section{
            .a1 = {-0x1.99cb3a3f35f6cp-2, -0x1.67211adbfd55ap-57},
            .a2 = {0x1.f9d4743f2ddf5p-2, -0x1.7fddec3773b9ap-61},
            .b0 = {0x1.18024e7ffdfa2p-2, 0x1.a038419822dbbp-57},
        },
        exact_section{
            .a1 = {-0x1.6c9db86fd1829p-2, -0x1.e7383d4f08c4ep-56},
            .a2 = {0x1.383bdc84a1aa0p-4, 0x1.cf157bfb1030ep-59},
            .b0 = {0x1.70b89f58ab73fp-3, 0x1.295539182ca0ap-57},
        },
    };
    constexpr double_double dc_target{0x1.fffe7db14b8a2p-1, 0x1.057f05701de2ap-58};

    auto const design = ctrlpp::make_chebyshev1<4>(100.0f, 1000.0f, 1e-4f);
    REQUIRE(design.has_value());
    check_against_exact_design(*design, exact, dc_target, counted_margins{.offset = 307.0, .numerator = 41.0, .dc_gain = 15.0});
}

TEST_CASE("make_chebyshev1 rejects a design whose coefficients degenerate",
          "[biquad][hardening][error]")
{
    // The ripple guard is a domain bound on the specification; the coefficient
    // sweep before the success return is what catches the remaining parameter
    // combinations that reach a non-finite coefficient through a path the
    // domain bounds do not describe. The case below satisfies every input
    // domain check and still degenerates, so it reaches the sweep, not a guard.
    SECTION("a positive ripple so small the ripple factor underflows to zero")
    {
        // ln(10) ripple_db / 10 underflows to zero at the smallest subnormal,
        // so the radicand is exactly zero, eps is zero, and the following
        // asinh(1 / eps) takes an infinite argument. The specification is
        // strictly positive, so the ripple guard passes it.
        auto const design
            = ctrlpp::make_chebyshev1<4>(100.0, 1000.0, std::numeric_limits<double>::denorm_min());
        REQUIRE(!design.has_value());
        CHECK(design.error() == ctrlpp::dsp_error::non_finite_input);
    }

}

TEST_CASE("make_chebyshev1 designs at a sample rate near the top of the double range",
          "[biquad][hardening][precision]")
{
    // The unit-scale formation forms tan(pi cutoff / fs) and never 2 fs tan(...), so it does not overflow here.
    // With an even order the DC gain sits at the ripple's trough 1 / sqrt(1 + eps^2), which puts the passband peak
    // at 0 dB. The DC gain matches it within the normalization's 14 roundings, this evaluation's 12 and the target's
    // own chain counted twice, 7: 33 eps.
    auto const design = ctrlpp::make_chebyshev1<4>(4.9e299, 1e300, 1.0);
    REQUIRE(design.has_value());

    double dc = 1.0;
    for(std::size_t k = 0; k < 2; ++k)
    {
        CAPTURE(k);
        auto const& c = design->section(k).coefficients();
        REQUIRE(std::isfinite(c.b0));
        REQUIRE(std::isfinite(c.b1));
        REQUIRE(std::isfinite(c.b2));
        REQUIRE(std::isfinite(c.a1));
        REQUIRE(std::isfinite(c.a2));
        dc *= (c.b0 + c.b1 + c.b2) / (1.0 + c.a1 + c.a2);
    }
    double const ripple = std::sqrt(std::expm1(std::numbers::ln10 / 10.0));
    double const trough = 1.0 / std::sqrt(1.0 + ripple * ripple);
    CHECK_THAT(dc, WithinRel(trough, 33.0 * eps));
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
