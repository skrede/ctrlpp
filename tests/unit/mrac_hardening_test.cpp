#include "hardening_helpers.h"
#include "ctrlpp/control/mrac.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

Vec1 vec1(double v)
{
    Vec1 r;
    r << v;
    return r;
}

auto make_ref_model() -> ctrlpp::siso_discrete_state_space<double, 1>
{
    ctrlpp::siso_discrete_state_space<double, 1> ref{};
    ref.A << 0.9;
    ref.B << 0.1;
    ref.C << 1.0;
    ref.D << 0.0;
    return ref;
}

auto make_siso_config(double gamma) -> ctrlpp::mrac_controller<double>::config_type
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << gamma;
    cfg.gamma_r << gamma;
    return cfg;
}

}

TEST_CASE("MRAC NaN state is rejected without touching the adaptive parameters",
          "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(0.5);
    ctrlpp::mrac_controller<double> ctrl(cfg);
    // Adapt once first, so the parameter matrices hold something other than
    // their initial values and "unchanged" is a real claim.
    REQUIRE(ctrl.evaluate(vec1(0.5), vec1(1.0)).has_value());

    const auto theta_x_before = ctrl.theta_x();
    const auto theta_r_before = ctrl.theta_r();
    const auto x_model_before = ctrl.x_model();
    const auto tracking_error_before = ctrl.tracking_error();

    auto rejected = ctrl.evaluate(ctrlpp::test::nan_vector<double, 1>(), vec1(1.0));

    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error() == ctrlpp::mrac_step_error::non_finite_state);
    // Exact, never a tolerance: a rejected cycle performs no arithmetic on the
    // carried state, so bitwise equality is the contract. Every piece is
    // asserted, not a representative one.
    CHECK(theta_x_before == ctrl.theta_x());
    CHECK(theta_r_before == ctrl.theta_r());
    CHECK(x_model_before == ctrl.x_model());
    CHECK(tracking_error_before == ctrl.tracking_error());
    CHECK(ctrl.health() == ctrlpp::mrac_health::ok);

    // The parameter matrices are the controller's memory: nothing re-derives
    // them, so had the poisoned sample been admitted, no later sample would
    // have recovered them. Compare against an instance that never saw it, and
    // compare the PARAMETERS, not only the command -- a command that happens to
    // match while the parameters diverged is the failure this catches.
    ctrlpp::mrac_controller<double> reference(cfg);
    REQUIRE(reference.evaluate(vec1(0.5), vec1(1.0)).has_value());

    auto after = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.25), vec1(1.0)));
    auto expected = ctrlpp::test::commanded(reference.evaluate(vec1(0.25), vec1(1.0)));
    CHECK(after == expected);
    CHECK(ctrl.theta_x() == reference.theta_x());
    CHECK(ctrl.theta_r() == reference.theta_r());
    CHECK(ctrl.x_model() == reference.x_model());
}

TEST_CASE("MRAC NaN reference is rejected without touching the adaptive parameters",
          "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(0.5);
    ctrlpp::mrac_controller<double> ctrl(cfg);
    REQUIRE(ctrl.evaluate(vec1(0.5), vec1(1.0)).has_value());

    const auto theta_x_before = ctrl.theta_x();
    const auto theta_r_before = ctrl.theta_r();
    const auto x_model_before = ctrl.x_model();
    const auto tracking_error_before = ctrl.tracking_error();

    auto rejected = ctrl.evaluate(vec1(0.0), ctrlpp::test::nan_vector<double, 1>());

    REQUIRE_FALSE(rejected.has_value());
    // A bad reference names the command generator, a bad state names the sensor
    // or estimator. Different subsystems, so different enumerators.
    CHECK(rejected.error() == ctrlpp::mrac_step_error::non_finite_reference);
    CHECK(theta_x_before == ctrl.theta_x());
    CHECK(theta_r_before == ctrl.theta_r());
    CHECK(x_model_before == ctrl.x_model());
    CHECK(tracking_error_before == ctrl.tracking_error());
    CHECK(ctrl.health() == ctrlpp::mrac_health::ok);

    ctrlpp::mrac_controller<double> reference(cfg);
    REQUIRE(reference.evaluate(vec1(0.5), vec1(1.0)).has_value());

    auto after = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.25), vec1(1.0)));
    auto expected = ctrlpp::test::commanded(reference.evaluate(vec1(0.25), vec1(1.0)));
    CHECK(after == expected);
    CHECK(ctrl.theta_x() == reference.theta_x());
    CHECK(ctrl.theta_r() == reference.theta_r());
    CHECK(ctrl.x_model() == reference.x_model());
}

TEST_CASE("MRAC zero gamma produces zero adaptation", "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(0.0);
    ctrlpp::mrac_controller<double> ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 100; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x), vec1(1.0)));
        x = 0.8 * x + 0.5 * u[0];
    }

    // With zero gamma, theta_x and theta_r should remain at initial (zero)
    REQUIRE_THAT(ctrl.theta_x()(0, 0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.theta_r()(0, 0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("MRAC huge gamma 1e15 produces finite output", "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(1e15);
    ctrlpp::mrac_controller<double> ctrl(cfg);

    auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));
    // First step with x=0: theta_x*0 + theta_r*1 = 0, but adaptation runs
    CHECK(std::isfinite(u[0]));
}

TEST_CASE("MRAC known first-order gain after one step", "[mrac][hardening][precision]")
{
    auto cfg = make_siso_config(0.5);
    ctrlpp::mrac_controller<double> ctrl(cfg);

    // Step 1: x=0, r=1
    // x_model advances: 0.9*0 + 0.1*1 = 0.1
    // tracking_error = x - x_model = 0 - 0.1 = -0.1
    // e_proj = B^T * tracking_error = 0.1 * (-0.1) = -0.01
    // theta_x -= sign_b * e_proj * x^T * gamma_x = 1 * (-0.01) * 0 * 0.5 = 0
    // theta_r -= sign_b * e_proj * r^T * gamma_r = 1 * (-0.01) * 1 * 0.5 = 0.005
    REQUIRE(ctrl.evaluate(vec1(0.0), vec1(1.0)).has_value());

    REQUIRE_THAT(ctrl.theta_x()(0, 0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.theta_r()(0, 0), WithinAbs(0.005, 1e-14));
}

TEST_CASE("MRAC with dead-zone tracks step reference within 5%", "[mrac][hardening][convergence]")
{
    using DeadZoneMrac = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::dead_zone>;
    DeadZoneMrac::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    cfg.robustification.threshold = 0.01;

    DeadZoneMrac ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 5000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x), vec1(1.0)));
        x = 0.8 * x + 0.5 * u[0];
    }

    auto x_m = ctrl.x_model()[0];
    REQUIRE(std::abs(x - x_m) < 0.05 * std::abs(x_m));
}

TEST_CASE("MRAC sigma-modification with moderate gamma keeps parameters bounded",
          "[mrac][hardening][robustness]")
{
    using SigmaMrac = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::sigma_modification>;
    SigmaMrac::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 10.0;
    cfg.gamma_r << 10.0;
    cfg.robustification.sigma = 1.0;

    SigmaMrac ctrl(cfg);

    double x = 0.0;
    bool all_bounded = true;

    for(int k = 0; k < 2000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x), vec1(1.0)));
        x = 0.8 * x + 0.5 * u[0];

        if(std::abs(ctrl.theta_x()(0, 0)) > 1e6 || std::abs(ctrl.theta_r()(0, 0)) > 1e6)
        {
            all_bounded = false;
            break;
        }
    }

    REQUIRE(all_bounded);
}

TEST_CASE("MRAC e-modification with moderate gamma keeps parameters bounded",
          "[mrac][hardening][robustness]")
{
    using EmodMrac = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::e_modification>;
    EmodMrac::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 10.0;
    cfg.gamma_r << 10.0;
    cfg.robustification.delta = 1.0;

    EmodMrac ctrl(cfg);

    double x = 0.0;
    bool all_bounded = true;

    for(int k = 0; k < 2000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x), vec1(1.0)));
        x = 0.8 * x + 0.5 * u[0];

        if(std::abs(ctrl.theta_x()(0, 0)) > 1e6 || std::abs(ctrl.theta_r()(0, 0)) > 1e6)
        {
            all_bounded = false;
            break;
        }
    }

    REQUIRE(all_bounded);
}
