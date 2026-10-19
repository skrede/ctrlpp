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
// Formulae for Audio EQ Biquad Filter Coefficients", 2005. With the
// Butterworth quality factor Q = 1/sqrt(2), `biquad::low_pass` meets both:
// -3.01 dB at the cutoff and a monotone passband.
//
// A "dirty derivative" (a differentiator low-pass filtered above its own
// bandwidth) behaves, at frequencies well below that bandwidth, like an
// ideal differentiator: |H(f)| ~ 2*pi*f. This is the textbook low-frequency
// limit of a first-order high-pass-shaped differentiator, not a property
// derived from this codebase. With the wc numerator factor in place,
// `biquad::dirty_derivative` tracks that 2*pi*f limit.

#include "ctrlpp/dsp/biquad.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <numbers>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

/// Exact digital magnitude |H(e^{jOmega})| at Omega = 2*pi*f/fs, evaluated
/// in double directly from the biquad's own reported coefficients.
template <typename Scalar>
auto magnitude_linear(const biquad_coeffs<Scalar>& c, double f, double fs) -> double
{
    const double omega = 2.0 * std::numbers::pi * f / fs;
    const std::complex<double> z_inv = std::exp(std::complex<double>(0.0, -omega));
    const std::complex<double> num = static_cast<double>(c.b0) + static_cast<double>(c.b1) * z_inv + static_cast<double>(c.b2) * z_inv * z_inv;
    const std::complex<double> den = 1.0 + static_cast<double>(c.a1) * z_inv + static_cast<double>(c.a2) * z_inv * z_inv;
    return std::abs(num / den);
}

/// |1 + a1 z^-1 + a2 z^-2| at Omega = 2*pi*f/fs: the denominator a coefficient perturbation moves.
auto denominator_magnitude(const biquad_coeffs<double>& c, double f, double fs) -> double
{
    const double omega = 2.0 * std::numbers::pi * f / fs;
    const std::complex<double> z_inv = std::exp(std::complex<double>(0.0, -omega));
    return std::abs(1.0 + c.a1 * z_inv + c.a2 * z_inv * z_inv);
}

auto magnitude_db(const biquad_coeffs<double>& c, double f, double fs) -> double
{
    return 20.0 * std::log10(magnitude_linear(c, f, fs));
}

// Running error bound (Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Sec. 3.3): a value
// carries how many eps its relative error can reach to first order, and every operation adds one rounding, exact or not.
struct counted
{
    double value;
    double roundings;
};

auto exact(double v) -> counted
{
    return {v, 0.0};
}

auto operator*(counted a, counted b) -> counted
{
    return {a.value * b.value, a.roundings + b.roundings + 1.0};
}

auto operator/(counted a, counted b) -> counted
{
    return {a.value / b.value, a.roundings + b.roundings + 1.0};
}

auto operator+(counted a, counted b) -> counted
{
    const double sum = a.value + b.value;
    return {sum, (std::abs(a.value) * a.roundings + std::abs(b.value) * b.roundings) / std::abs(sum) + 1.0};
}

/// A library function with a one-ulp error; slope = |x f'(x) / f(x)| carries the argument's relative error through it.
auto through(double fx, double slope, counted x) -> counted
{
    return {fx, slope * x.roundings + 1.0};
}

/// Omega = tan(pi * (fc / fs)) as biquad.h forms it, pi rounded into the scalar type.
auto count_prewarp(double fc, double fs) -> counted
{
    const counted x = counted{std::numbers::pi, 1.0} * (exact(fc) / exact(fs));
    return through(std::tan(x.value), 2.0 * x.value / std::sin(2.0 * x.value), x);
}

struct section_count
{
    counted b0;
    counted q1;
    counted p2;
    counted dc_denominator;
};

/// The offset-form section as biquad.h forms it: b0, the a1 offset q1 = a1 + 2, the a2 offset p2 = 1 - a2, and D.
auto count_section(counted omega, counted g1, counted g0) -> section_count
{
    const counted beta1 = omega * g1;
    const counted beta0 = omega * omega * g0;
    const counted den0 = exact(1.0) + beta1 + beta0;
    return {beta0 / den0, (exact(2.0) * beta1 + exact(4.0) * beta0) / den0, exact(2.0) * beta1 / den0, exact(4.0) * beta0 / den0};
}

auto count_low_pass(double fc, double fs) -> section_count
{
    return count_section(count_prewarp(fc, fs), counted{std::numbers::sqrt2, 1.0}, exact(1.0));
}

constexpr double eps_single = static_cast<double>(std::numeric_limits<float>::epsilon());
constexpr double eps_double = std::numeric_limits<double>::epsilon();

// Storing a1 in [1, 2) costs half its ulp, eps / 2, and a2 below 1 at most eps / 4.
constexpr double offset_storage_roundings = 0.75;

auto float_ulp(float v) -> double
{
    const float a = std::abs(v);
    return static_cast<double>(std::nextafter(a, std::numeric_limits<float>::infinity()) - a);
}

/// |float coefficient - double design rounded to float|: half an ulp for each of the two stores into float, the
/// counted formation error of the term the coefficient is stored from at both precisions, and the double store.
void check_coefficient(float single, double reference, counted term)
{
    const float rounded = static_cast<float>(reference);
    const double error = std::abs(static_cast<double>(single) - static_cast<double>(rounded));
    const double bound = float_ulp(rounded) + (eps_single + eps_double) * term.roundings * std::abs(term.value) + eps_double * std::abs(reference);
    CAPTURE(single, reference, error / float_ulp(rounded), bound / float_ulp(rounded));
    CHECK(error <= bound);
}

void check_section_coefficients(const biquad_coeffs<float>& single, const biquad_coeffs<double>& reference, const section_count& n)
{
    check_coefficient(single.b0, reference.b0, n.b0);
    check_coefficient(single.b1, reference.b1, counted{2.0 * n.b0.value, n.b0.roundings});
    check_coefficient(single.b2, reference.b2, n.b0);
    check_coefficient(single.a1, reference.a1, n.q1);
    check_coefficient(single.a2, reference.a2, n.p2);
}

/// Relative error of a stored float section's response where its denominator 1 + a1 z^-1 + a2 z^-2 has magnitude
/// denominator: the numerator's counted error, and the denominator moved by at most the offset storage plus the
/// counted formation of the two offsets. Exact in the denominator perturbation, which nears 1 at the floor.
auto section_response_bound(const section_count& n, double denominator) -> double
{
    const double moved = (offset_storage_roundings + n.q1.roundings * n.q1.value + n.p2.roundings * n.p2.value) * eps_single / denominator;
    REQUIRE(moved < 1.0);
    return (n.b0.roundings * eps_single + moved) / (1.0 - moved);
}

// Evaluating the response in double from the coefficients: the rounding-op margin of the anchors above.
constexpr double evaluation_bound = 16.0 * eps_double;

/// Butterworth section k of an order-Order cascade: g1 = 2 cos(theta_k), theta_k = pi (2k + 1) / (2 Order).
auto count_butterworth_section(counted omega, std::size_t k, std::size_t order) -> section_count
{
    const counted theta = counted{std::numbers::pi, 1.0} * exact(static_cast<double>(2 * k + 1)) / exact(static_cast<double>(2 * order));
    const counted cosine = through(std::cos(theta.value), theta.value * std::tan(theta.value), theta);
    return count_section(omega, exact(2.0) * cosine, exact(1.0));
}

struct chebyshev1_count
{
    std::array<section_count, 2> sections;
    counted target_dc;
};

/// make_chebyshev1<4>'s chain: the ripple factor through expm1 and sqrt, the prototype through asinh, sinh and cosh,
/// each section's pole pair through sin and cos, and the DC target 1 / sqrt(1 + eps^2).
auto count_chebyshev1_order4(double fc, double fs, double ripple_db) -> chebyshev1_count
{
    const counted y = counted{std::numbers::ln10, 1.0} * exact(ripple_db) / exact(10.0);
    const counted radicand = through(std::expm1(y.value), y.value * std::exp(y.value) / std::expm1(y.value), y);
    const counted ripple = through(std::sqrt(radicand.value), 0.5, radicand);
    const counted z = exact(1.0) / ripple;
    const double arc = std::asinh(z.value);
    const counted v = through(arc, z.value / (std::sqrt(1.0 + z.value * z.value) * arc), z) / exact(4.0);
    const counted sinh_v = through(std::sinh(v.value), v.value / std::tanh(v.value), v);
    const counted cosh_v = through(std::cosh(v.value), v.value * std::tanh(v.value), v);
    const counted omega = count_prewarp(fc, fs);

    chebyshev1_count n{};
    for(std::size_t k = 0; k < 2; ++k)
    {
        const counted theta = counted{std::numbers::pi, 1.0} * exact(static_cast<double>(2 * k + 1)) / exact(8.0);
        const counted sigma = sinh_v * through(std::sin(theta.value), theta.value / std::tan(theta.value), theta);
        const counted imag = cosh_v * through(std::cos(theta.value), theta.value * std::tan(theta.value), theta);
        n.sections[k] = count_section(omega, exact(2.0) * sigma, sigma * sigma + imag * imag);
    }
    const counted root = through(std::sqrt(1.0 + ripple.value * ripple.value), 0.5, exact(1.0) + ripple * ripple);
    n.target_dc = exact(1.0) / root;
    return n;
}

// The DC normalization of a two-section cascade: per section two numerator and two denominator sums, the quotient
// and the running product, then the correction quotient and the scaling of the first section's numerator.
constexpr double chebyshev1_normalization_roundings = 2.0 * 6.0 + 2.0;

/// Bounds on |a1_float - a1_double| and |a2_float - a2_double|: the half-ulp store of a1 below 2 and of a2 below 1 at
/// either precision, and the counted formation of both offsets at either precision.
auto offset_perturbation(const section_count& n) -> std::array<double, 2>
{
    const double both = eps_single + eps_double;
    return {both / 2.0 + n.q1.roundings * n.q1.value * both, both / 4.0 + n.p2.roundings * n.p2.value * both};
}

/// | |H_f(fc) / H_f(0)| / |H_d(fc) / H_d(0)| - 1 | for a section whose numerator shape is fixed and whose a1, a2 move by
/// at most s[0], s[1]. With w = (d1 z^-1 + d2 z^-2) / A(z) the ratio is (1 + (d1 + d2) / gap) / |1 + w|, and
/// 1 + Re w <= |1 + w| <= 1 + Re w + (Im w)^2 / (2 (1 + Re w)) bounds it to all orders.
auto shape_ratio_bound(const biquad_coeffs<double>& d, double gap, double f, double fs, std::array<double, 2> s) -> double
{
    const std::complex<double> z_inv = std::exp(std::complex<double>(0.0, -2.0 * std::numbers::pi * f / fs));
    const std::complex<double> a = 1.0 + d.a1 * z_inv + d.a2 * z_inv * z_inv;
    const std::complex<double> w1 = z_inv / a;
    const std::complex<double> w2 = z_inv * z_inv / a;
    const double linear = s[0] * std::abs(1.0 / gap - w1.real()) + s[1] * std::abs(1.0 / gap - w2.real());
    const double real = s[0] * std::abs(w1.real()) + s[1] * std::abs(w2.real());
    const double imag = s[0] * std::abs(w1.imag()) + s[1] * std::abs(w2.imag());
    REQUIRE(real < 1.0);
    return (linear + imag * imag / (2.0 * (1.0 - real))) / (1.0 - real);
}

auto dc_gain(const biquad_coeffs<float>& c) -> double
{
    const double numerator = static_cast<double>(c.b0) + static_cast<double>(c.b1) + static_cast<double>(c.b2);
    return numerator / (1.0 + static_cast<double>(c.a1) + static_cast<double>(c.a2));
}

/// Brackets the floor of a unit-scale second-order section with g0 = 1: the design closed over its frequency ratio is
/// accepted at r* (1 + m eps) and refused at r* (1 - m eps), where D = count * eps.
template <typename Scalar, typename Design>
void check_second_order_floor_boundary(long double count, counted g1, Design design)
{
    using wide = long double;
    const wide eps = static_cast<wide>(std::numeric_limits<Scalar>::epsilon());
    const wide c_eps = count * eps;
    // 4 Omega^2 / (1 + g1 Omega + Omega^2) = c eps solved for its positive root, no cancellation in the numerator.
    const wide a = wide{4} - c_eps;
    const wide b = c_eps * static_cast<wide>(g1.value);
    const wide omega = (b + std::sqrt(b * b + wide{4} * a * c_eps)) / (wide{2} * a);
    const wide ratio = std::atan(omega) / std::numbers::pi_v<wide>;

    // The D chain at Scalar precision plus rounding the ratio into Scalar, and the closed form's at most fourteen
    // roundings above at long double precision expressed in Scalar eps. D grows as the ratio squared, so a ratio m eps
    // off r* moves D twice as far as its own computed error.
    const wide d_chain = static_cast<wide>(count_section(count_prewarp(static_cast<double>(ratio), 1.0), g1, exact(1.0)).dc_denominator.roundings + 1.0);
    const wide m = d_chain + wide{14} * std::numeric_limits<wide>::epsilon() / eps;
    CAPTURE(static_cast<double>(ratio), static_cast<double>(m));

    const auto above = design(static_cast<Scalar>(ratio * (wide{1} + m * eps)));
    REQUIRE(above.has_value());
    const auto below = design(static_cast<Scalar>(ratio * (wide{1} - m * eps)));
    REQUIRE(!below.has_value());
    CHECK(below.error() == dsp_error::cutoff_below_coefficient_resolution);
}

template <typename Scalar>
void check_low_pass_floor_boundary()
{
    check_second_order_floor_boundary<Scalar>(3.0L / 4.0L, counted{std::numbers::sqrt2, 1.0},
                                              [](Scalar ratio) { return biquad<Scalar>::low_pass(ratio, Scalar{1}); });
}

template <typename Scalar>
void check_notch_floor_boundary(Scalar q)
{
    check_second_order_floor_boundary<Scalar>(7.0L / 4.0L, exact(1.0) / exact(static_cast<double>(q)),
                                              [q](Scalar ratio) { return biquad<Scalar>::notch(ratio, Scalar{1}, q); });
}

auto count_notch(double f, double fs, double q) -> section_count
{
    return count_section(count_prewarp(f, fs), exact(1.0) / exact(q), exact(1.0));
}

struct notch_errors
{
    double b0;
    double a1;
    double a2;
};

/// Absolute errors, in eps of the precision a notch is formed in, of b0 = b2, b1 = a1 and a2: the store of each (b0
/// below 1 at eps / 4, a1 below 2 at eps / 2, a2 below 1 at eps / 4) and the counted formation of the offset it is
/// stored from, b0 = 1 - p2 / 2, a1 = -2 + q1, a2 = 1 - p2.
auto notch_coefficient_errors(const section_count& n) -> notch_errors
{
    return {0.25 + n.p2.roundings * n.p2.value / 2.0, 0.5 + n.q1.roundings * n.q1.value, 0.25 + n.p2.roundings * n.p2.value};
}

void check_notch_coefficients(const biquad_coeffs<float>& single, const biquad_coeffs<double>& reference, const section_count& n)
{
    const counted half_offset{n.p2.value / 2.0, n.p2.roundings};
    check_coefficient(single.b0, reference.b0, half_offset);
    check_coefficient(single.b1, reference.b1, n.q1);
    check_coefficient(single.b2, reference.b2, half_offset);
    check_coefficient(single.a1, reference.a1, n.q1);
    check_coefficient(single.a2, reference.a2, n.p2);
}

/// |H_stored / H - 1| for a numerator and a denominator moved by at most the relative errors e_n and e_a: exact in the
/// denominator perturbation, which nears 1 at the floor.
auto response_bound(double e_n, double e_a) -> double
{
    REQUIRE(e_a < 1.0);
    return (e_n + e_a) / (1.0 - e_a);
}

/// The notch response bound at a frequency where the design's numerator and denominator have the given magnitudes.
auto notch_response_bound(const section_count& n, double numerator, double denominator, double eps) -> double
{
    const notch_errors e = notch_coefficient_errors(n);
    return response_bound((2.0 * e.b0 + e.a1) * eps / numerator, (e.a1 + e.a2) * eps / denominator);
}

auto numerator_magnitude(const biquad_coeffs<double>& c, double f, double fs) -> double
{
    const std::complex<double> z_inv = std::exp(std::complex<double>(0.0, -2.0 * std::numbers::pi * f / fs));
    return std::abs(c.b0 + c.b1 * z_inv + c.b2 * z_inv * z_inv);
}

struct differentiator_count
{
    counted b0;
    counted gap;
};

/// The first-order section as biquad.h forms it: share = Omega / (1 + Omega), b0 = 2 fs share, D1 = 2 share.
auto count_differentiator(double f, double fs) -> differentiator_count
{
    const counted omega = count_prewarp(f, fs);
    const counted share = omega / (exact(1.0) + omega);
    return {exact(2.0 * fs) * share, exact(2.0) * share};
}

/// Relative errors of a dirty derivative's numerator b0 (1 - z^-1) and of its a1 = -1 + D1 stored below 1 at eps / 4.
auto differentiator_response_bound(const differentiator_count& n, double denominator, double eps) -> double
{
    return response_bound(n.b0.roundings * eps, (0.25 + n.gap.roundings * n.gap.value) * eps / denominator);
}

template <typename Scalar>
void check_first_order_floor_boundary()
{
    using wide = long double;
    const wide eps = static_cast<wide>(std::numeric_limits<Scalar>::epsilon());
    const wide c_eps = eps / wide{4};
    // 2 Omega / (1 + Omega) = c eps solved for Omega.
    const wide omega = c_eps / (wide{2} - c_eps);
    const wide ratio = std::atan(omega) / std::numbers::pi_v<wide>;

    // The D1 chain at Scalar precision plus rounding the ratio into Scalar, and the closed form's five roundings above
    // (two for Omega, atan, pi and the division by it). D1 grows as the ratio, so m is the computed error itself.
    const wide d_chain = static_cast<wide>(count_differentiator(static_cast<double>(ratio), 1.0).gap.roundings + 1.0);
    const wide m = d_chain + wide{5} * std::numeric_limits<wide>::epsilon() / eps;
    CAPTURE(static_cast<double>(ratio), static_cast<double>(m));

    const auto above = biquad<Scalar>::dirty_derivative(static_cast<Scalar>(ratio * (wide{1} + m * eps)), Scalar{1});
    REQUIRE(above.has_value());
    const auto below = biquad<Scalar>::dirty_derivative(static_cast<Scalar>(ratio * (wide{1} - m * eps)), Scalar{1});
    REQUIRE(!below.has_value());
    CHECK(below.error() == dsp_error::cutoff_below_coefficient_resolution);
}

} // namespace

TEST_CASE("Butterworth low-pass biquad reaches -3.01 dB at cutoff with a monotone passband", "[dsp][anchor]")
{
    const double fc = 100.0;
    const double fs = 1000.0;
    const auto lp = biquad<double>::low_pass(fc, fs);
    REQUIRE(lp.has_value());
    const auto c = lp->coefficients();

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

TEST_CASE("dirty_derivative magnitude tracks 2*pi*f well below its bandwidth", "[dsp][anchor]")
{
    const double bandwidth_hz = 50.0;
    const double fs = 1000.0;
    const auto dd = biquad<double>::dirty_derivative(bandwidth_hz, fs);
    REQUIRE(dd.has_value());
    const auto c = dd->coefficients();

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

TEST_CASE("Single-precision low-pass designs carry the double design down to the coefficient-resolution floor", "[dsp][anchor][float]")
{
    constexpr std::array<double, 11> accepted{1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 5e-4, 3e-4, 2e-4, 1e-4, 7e-5, 5e-5};
    constexpr float fs = 1000.0f;

    for(double ratio : accepted)
    {
        CAPTURE(ratio);
        const float fc = static_cast<float>(ratio * static_cast<double>(fs));
        const auto single = biquad<float>::low_pass(fc, fs);
        const auto reference = biquad<double>::low_pass(static_cast<double>(fc), static_cast<double>(fs));
        REQUIRE(single.has_value());
        REQUIRE(reference.has_value());
        const auto n = count_low_pass(static_cast<double>(fc), static_cast<double>(fs));
        check_section_coefficients(single->coefficients(), reference->coefficients(), n);

        const double dc_bound = section_response_bound(n, n.dc_denominator.value) + evaluation_bound;
        CHECK(std::abs(dc_gain(single->coefficients()) - 1.0) <= dc_bound);

        const double at_cutoff = magnitude_linear(single->coefficients(), static_cast<double>(fc), static_cast<double>(fs));
        const double denominator = denominator_magnitude(reference->coefficients(), static_cast<double>(fc), static_cast<double>(fs));
        const double cutoff_bound = section_response_bound(n, denominator) + evaluation_bound;
        CAPTURE(at_cutoff, dc_bound, cutoff_bound);
        CHECK(std::abs(at_cutoff * std::numbers::sqrt2 - 1.0) <= cutoff_bound);
    }

    for(double ratio : {3e-5, 1e-5})
    {
        CAPTURE(ratio);
        const auto single = biquad<float>::low_pass(static_cast<float>(ratio * static_cast<double>(fs)), fs);
        REQUIRE(!single.has_value());
        CHECK(single.error() == dsp_error::cutoff_below_coefficient_resolution);
    }
}

TEST_CASE("The low-pass coefficient-resolution floor sits where its closed form puts it", "[dsp][anchor][float]")
{
    check_low_pass_floor_boundary<float>();
    check_low_pass_floor_boundary<double>();
}

TEST_CASE("A single-precision Butterworth cascade at 20 Hz / 48 kHz carries the double design", "[dsp][anchor][float]")
{
    constexpr float fc = 20.0f;
    constexpr float fs = 48000.0f;
    const auto single = make_butterworth<4>(fc, fs);
    const auto reference = make_butterworth<4>(static_cast<double>(fc), static_cast<double>(fs));
    REQUIRE(single.has_value());
    REQUIRE(reference.has_value());

    const counted omega = count_prewarp(static_cast<double>(fc), static_cast<double>(fs));
    double at_cutoff = 1.0;
    double dc = 1.0;
    double cutoff_bound = 1.0;
    double dc_bound = 1.0;
    for(std::size_t k = 0; k < 2; ++k)
    {
        CAPTURE(k);
        const auto n = count_butterworth_section(omega, k, 4);
        const auto& c = single->section(k).coefficients();
        const auto& d = reference->section(k).coefficients();
        check_section_coefficients(c, d, n);
        at_cutoff *= magnitude_linear(c, static_cast<double>(fc), static_cast<double>(fs));
        dc *= dc_gain(c);
        cutoff_bound *= 1.0 + section_response_bound(n, denominator_magnitude(d, static_cast<double>(fc), static_cast<double>(fs))) + evaluation_bound;
        dc_bound *= 1.0 + section_response_bound(n, n.dc_denominator.value) + evaluation_bound;
    }
    CAPTURE(20.0 * std::log10(at_cutoff), dc, cutoff_bound - 1.0, dc_bound - 1.0);
    CHECK(std::abs(at_cutoff * std::numbers::sqrt2 - 1.0) <= cutoff_bound - 1.0);
    CHECK(std::abs(dc - 1.0) <= dc_bound - 1.0);
}

TEST_CASE("A single-precision Chebyshev cascade at 20 Hz / 48 kHz matches the double design within its representation bound",
          "[dsp][anchor][float]")
{
    constexpr float fc = 20.0f;
    constexpr float fs = 48000.0f;
    const auto single = make_chebyshev1<4>(fc, fs, 1.0f);
    const auto reference = make_chebyshev1<4>(static_cast<double>(fc), static_cast<double>(fs), 1.0);
    REQUIRE(single.has_value());
    REQUIRE(reference.has_value());
    const auto n = count_chebyshev1_order4(static_cast<double>(fc), static_cast<double>(fs), 1.0);

    // The normalization fixes each cascade's DC gain to its own target, and every section's numerator shape is
    // (1 + z^-1)^2 in both, so the two responses at the cutoff differ only through the stored denominators.
    const double normalization = (n.target_dc.roundings + chebyshev1_normalization_roundings) * (eps_single + eps_double);
    double upper = 1.0 + normalization;
    double lower = 1.0 - normalization;
    double single_cutoff = 1.0;
    double reference_cutoff = 1.0;
    for(std::size_t k = 0; k < 2; ++k)
    {
        CAPTURE(k);
        const auto& c = single->section(k).coefficients();
        const auto& d = reference->section(k).coefficients();
        const double gap = n.sections[k].dc_denominator.value;
        const double shape = shape_ratio_bound(d, gap, static_cast<double>(fc), static_cast<double>(fs), offset_perturbation(n.sections[k]));
        upper *= 1.0 + shape;
        lower *= 1.0 - shape;
        // The first section's numerator carries the normalization, which absorbs the cascade's DC error.
        if(k == 0)
        {
            check_coefficient(c.a1, d.a1, n.sections[k].q1);
            check_coefficient(c.a2, d.a2, n.sections[k].p2);
        }
        else
        {
            check_section_coefficients(c, d, n.sections[k]);
        }
        single_cutoff *= magnitude_linear(c, static_cast<double>(fc), static_cast<double>(fs));
        reference_cutoff *= magnitude_linear(d, static_cast<double>(fc), static_cast<double>(fs));
    }
    const double ratio = single_cutoff / reference_cutoff;
    const double evaluation = 4.0 * evaluation_bound;
    CAPTURE(20.0 * std::log10(single_cutoff), 20.0 * std::log10(reference_cutoff), ratio - 1.0, upper - 1.0, 1.0 - lower);
    CHECK(ratio <= upper + evaluation);
    CHECK(ratio >= lower - evaluation);
}

TEST_CASE("Single-precision notch designs carry the double design down to the coefficient-resolution floor", "[dsp][anchor][float]")
{
    constexpr std::array<double, 5> accepted{1e-1, 1e-2, 1e-3, 3e-4, 1e-4};
    constexpr float fs = 1000.0f;
    constexpr float q = 0.7f;
    const double quarter = static_cast<double>(fs) / 4.0;

    for(double ratio : accepted)
    {
        CAPTURE(ratio);
        const float f = static_cast<float>(ratio * static_cast<double>(fs));
        const auto single = biquad<float>::notch(f, fs, q);
        const auto reference = biquad<double>::notch(static_cast<double>(f), static_cast<double>(fs), static_cast<double>(q));
        REQUIRE(single.has_value());
        REQUIRE(reference.has_value());
        const auto n = count_notch(static_cast<double>(f), static_cast<double>(fs), static_cast<double>(q));
        check_notch_coefficients(single->coefficients(), reference->coefficients(), n);

        // The numerator gap equals the design gap D, so the DC gain is 1.
        const double gap = n.dc_denominator.value;
        const double dc_bound = notch_response_bound(n, gap, gap, eps_single) + evaluation_bound;
        CHECK(std::abs(dc_gain(single->coefficients()) - 1.0) <= dc_bound);

        const double numerator = numerator_magnitude(reference->coefficients(), quarter, static_cast<double>(fs));
        const double denominator = denominator_magnitude(reference->coefficients(), quarter, static_cast<double>(fs));
        const double e_single = notch_response_bound(n, numerator, denominator, eps_single);
        const double e_double = notch_response_bound(n, numerator, denominator, eps_double);
        const double quarter_bound = (e_single + e_double) / (1.0 - e_double) + 2.0 * evaluation_bound;
        const double quarter_ratio = magnitude_linear(single->coefficients(), quarter, static_cast<double>(fs))
                                     / magnitude_linear(reference->coefficients(), quarter, static_cast<double>(fs));
        CAPTURE(dc_gain(single->coefficients()), dc_bound, quarter_ratio, quarter_bound);
        CHECK(std::abs(quarter_ratio - 1.0) <= quarter_bound);
    }

    for(double ratio : {7e-5, 3e-5, 1e-5})
    {
        CAPTURE(ratio);
        const auto single = biquad<float>::notch(static_cast<float>(ratio * static_cast<double>(fs)), fs, q);
        REQUIRE(!single.has_value());
        CHECK(single.error() == dsp_error::cutoff_below_coefficient_resolution);
    }
}

TEST_CASE("The notch coefficient-resolution floor sits where its closed form puts it", "[dsp][anchor][float]")
{
    check_notch_floor_boundary<float>(0.7f);
    check_notch_floor_boundary<double>(0.7);
}

TEST_CASE("A single-precision dirty derivative carries the double design down to its first-order floor", "[dsp][anchor][float]")
{
    constexpr std::array<double, 7> accepted{1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    constexpr float fs = 1000.0f;

    for(double ratio : accepted)
    {
        CAPTURE(ratio);
        const float f = static_cast<float>(ratio * static_cast<double>(fs));
        const auto single = biquad<float>::dirty_derivative(f, fs);
        const auto reference = biquad<double>::dirty_derivative(static_cast<double>(f), static_cast<double>(fs));
        REQUIRE(single.has_value());
        REQUIRE(reference.has_value());
        const auto n = count_differentiator(static_cast<double>(f), static_cast<double>(fs));
        const auto& c = single->coefficients();
        const auto& d = reference->coefficients();
        check_coefficient(c.b0, d.b0, n.b0);
        check_coefficient(c.b1, d.b1, n.b0);
        check_coefficient(c.a1, d.a1, n.gap);
        CHECK(c.b2 == 0.0f);
        CHECK(c.a2 == 0.0f);

        const double denominator = denominator_magnitude(d, static_cast<double>(f), static_cast<double>(fs));
        const double e_single = differentiator_response_bound(n, denominator, eps_single);
        const double e_double = differentiator_response_bound(n, denominator, eps_double);
        const double bound = (e_single + e_double) / (1.0 - e_double) + 2.0 * evaluation_bound;
        const double at_bandwidth = magnitude_linear(c, static_cast<double>(f), static_cast<double>(fs))
                                    / magnitude_linear(d, static_cast<double>(f), static_cast<double>(fs));
        CAPTURE(at_bandwidth, bound);
        CHECK(std::abs(at_bandwidth - 1.0) <= bound);
    }
}

TEST_CASE("The dirty-derivative first-order floor sits where its closed form puts it", "[dsp][anchor][float]")
{
    check_first_order_floor_boundary<float>();
    check_first_order_floor_boundary<double>();
}
