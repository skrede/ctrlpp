// This regression pins two properties of a default-constructed L1
// configuration and controller. First, the projection bounds default to an
// unbounded range, so a default config adapts freely instead of clamping the
// uncertainty estimate to exactly zero (which silently disabled adaptation).
// Second, the controller's default output filter C(s) is the corrected
// maximally-flat (Butterworth) low-pass inherited from `biquad::low_pass`:
// -3.01 dB at the cutoff, not the +3 dB peaking of a Q = sqrt(2) section.

#include "hardening_helpers.h"
#include "ctrlpp/control/l1.h"
#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/vector_biquad.h"
#include "ctrlpp/types.h"

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

using Vec1 = Vector<double, 1>;

auto vec1(double v) -> Vec1
{
    Vec1 r;
    r << v;
    return r;
}

// The predictor pole (A_m = 0.9) deliberately mismatches the plant pole
// (0.8), so the prediction error drives the uncertainty estimate away from
// zero. A config that never clamps therefore adapts; the old zero-default
// bounds pinned the estimate at exactly zero. theta_min / theta_max are left
// at their config defaults.
auto make_adapting_config() -> l1_config<double, 1, 1>
{
    l1_config<double, 1, 1> cfg{};
    cfg.predictor_model.A << 0.9;
    cfg.predictor_model.B << 0.5;
    cfg.predictor_model.C << 1.0;
    cfg.predictor_model.D << 0.0;
    cfg.gamma << 10.0;
    return cfg;
}

auto run_sigma_hat(const l1_config<double, 1, 1>& cfg) -> double
{
    // create() is the only construction path and it is fallible; this config is
    // valid, so a rejection here is a test failure rather than a skip.
    auto created = l1_controller<double>::create(cfg, 15.0, 100.0);
    REQUIRE(created.has_value());
    auto& ctrl = *created;
    double x_plant = 0.0;
    for(int k = 0; k < 200; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(1.0)));
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }
    return ctrl.sigma_hat()[0];
}

// Exact digital magnitude |H(e^{jOmega})| at Omega = 2*pi*f/fs, evaluated
// directly from the biquad's own reported coefficients (Oppenheim &amp;
// Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Ch. 5).
auto magnitude_linear(const biquad_coeffs<double>& c, double f, double fs) -> double
{
    const double omega = 2.0 * std::numbers::pi * f / fs;
    const std::complex<double> z_inv = std::exp(std::complex<double>(0.0, -omega));
    const std::complex<double> num = c.b0 + c.b1 * z_inv + c.b2 * z_inv * z_inv;
    const std::complex<double> den = 1.0 + c.a1 * z_inv + c.a2 * z_inv * z_inv;
    return std::abs(num / den);
}

} // namespace

TEST_CASE("l1 default projection bounds are unbounded", "[l1]")
{
    l1_config<double, 1, 1> cfg{};
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE(cfg.theta_min[0] == -inf);
    REQUIRE(cfg.theta_max[0] == inf);
}

TEST_CASE("l1 default config adapts instead of clamping the estimate to zero", "[l1]")
{
    // Old default: theta_min = theta_max = 0 clamped the estimate every step,
    // pinning it at exactly zero (adaptation silently disabled).
    auto zero_bounds = make_adapting_config();
    zero_bounds.theta_min << 0.0;
    zero_bounds.theta_max << 0.0;
    REQUIRE(run_sigma_hat(zero_bounds) == 0.0);

    // The unbounded default adapts: the estimate leaves zero.
    const double sigma_default = run_sigma_hat(make_adapting_config());
    REQUIRE(sigma_default != 0.0);

    // Finite bounds set wide enough never to bind must reproduce the unbounded
    // default exactly. In both runs the projection evaluates max(x, lo) and
    // min(x, hi) with a bound that never binds, so it returns x unchanged: the
    // two runs execute bit-identical arithmetic, not merely close values.
    auto wide_bounds = make_adapting_config();
    wide_bounds.theta_min << -1.0e3;
    wide_bounds.theta_max << 1.0e3;
    REQUIRE(sigma_default == run_sigma_hat(wide_bounds));
}

TEST_CASE("l1 default output filter inherits the Butterworth low-pass", "[l1]")
{
    const double fc = 100.0;
    const double fs = 1000.0;

    // The default L1 output filter is vector_biquad<Scalar, NU>, whose low_pass
    // fills every channel from the scalar biquad::low_pass prototype. Confirm
    // the vector-valued default filter reproduces that prototype exactly: same
    // coefficients and same state evolution give bit-identical output.
    auto vfilter = vector_biquad<double, 1>::low_pass(fc, fs);
    auto sfilter = biquad<double>::low_pass(fc, fs);
    REQUIRE(vfilter.has_value());
    REQUIRE(sfilter.has_value());
    for(int k = 0; k < 64; ++k)
    {
        const double x = (k == 0) ? 1.0 : 0.5;
        const double yv = vfilter->process(vec1(x))[0];
        const double ys = sfilter->process(x);
        CAPTURE(k, yv, ys);
        REQUIRE(yv == ys);
    }

    // That prototype is a maximally-flat low-pass: |H(fc)| = 1/sqrt(2), i.e.
    // -3.01 dB at the half-power point. The rounding-op margin converts the
    // handful of complex multiply-adds and one division (each up to one ULP at
    // the O(1) magnitude scale) into a decibel tolerance via the local
    // derivative of 20*log10(x)/ln(10) at the target magnitude.
    const auto prototype = biquad<double>::low_pass(fc, fs);
    REQUIRE(prototype.has_value());
    const auto c = prototype->coefficients();
    constexpr double rounding_op_margin = 16.0;
    const double eps = std::numeric_limits<double>::epsilon();
    const double target_mag = 1.0 / std::numbers::sqrt2;
    const double db_per_unit_mag = 20.0 / (std::log(10.0) * target_mag);
    const double tol_db = db_per_unit_mag * rounding_op_margin * eps;

    // -10*log10(2) is the exact value of 20*log10(1/sqrt(2)), derived here
    // rather than typed as a rounded literal.
    const double target_db = -10.0 * std::log10(2.0);
    const double measured_db = 20.0 * std::log10(magnitude_linear(c, fc, fs));
    CAPTURE(measured_db, target_db);
    REQUIRE_THAT(measured_db, WithinAbs(target_db, tol_db));
}
