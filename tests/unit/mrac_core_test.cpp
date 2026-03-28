#include "ctrlpp/control/mrac.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

auto make_ref_model() -> ctrlpp::siso_discrete_state_space<double, 1>
{
    ctrlpp::siso_discrete_state_space<double, 1> ref{};
    ref.A << 0.9;
    ref.B << 0.1;
    ref.C << 1.0;
    ref.D << 0.0;
    return ref;
}

using Vec1 = ctrlpp::Vector<double, 1>;

Vec1 vec1(double v)
{
    Vec1 r;
    r << v;
    return r;
}

}

TEST_CASE("MRAC default construction and accessors", "[mrac]")
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma = 0.5;
    cfg.sign_b = 1.0;

    ctrlpp::mrac_controller<double> ctrl(cfg);

    REQUIRE_THAT(ctrl.theta_x(), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.theta_r(), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.x_model()[0], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.tracking_error()[0], WithinAbs(0.0, 1e-15));
}

TEST_CASE("MRAC tracks step reference within 5% of reference model", "[mrac]")
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma = 0.5;
    cfg.sign_b = 1.0;

    ctrlpp::mrac_controller<double> ctrl(cfg);

    double x = 0.0;

    for(int k = 0; k < 500; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u;
    }

    auto x_m = ctrl.x_model()[0];
    REQUIRE(std::abs(x - x_m) < 0.05 * std::abs(x_m));
}

TEST_CASE("MRAC evaluate returns scalar control signal", "[mrac]")
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma = 0.5;
    cfg.sign_b = 1.0;

    ctrlpp::mrac_controller<double> ctrl(cfg);

    auto u = ctrl.evaluate(vec1(0.0), vec1(1.0));
    REQUIRE(std::isfinite(u));
}

TEST_CASE("MRAC adaptation modifies theta_x and theta_r", "[mrac]")
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma = 0.5;
    cfg.sign_b = 1.0;

    ctrlpp::mrac_controller<double> ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 10; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u;
    }

    CHECK(ctrl.theta_x() != 0.0);
    CHECK(ctrl.theta_r() != 0.0);
}

TEST_CASE("MRAC reset restores initial state", "[mrac]")
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma = 0.5;
    cfg.sign_b = 1.0;
    cfg.theta_x_0 = 0.0;
    cfg.theta_r_0 = 0.0;

    ctrlpp::mrac_controller<double> ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 50; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u;
    }

    REQUIRE(ctrl.theta_x() != 0.0);

    ctrl.reset();

    REQUIRE_THAT(ctrl.theta_x(), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.theta_r(), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.x_model()[0], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.tracking_error()[0], WithinAbs(0.0, 1e-15));
}

TEST_CASE("MRAC reference model state advances each step", "[mrac]")
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma = 0.5;
    cfg.sign_b = 1.0;

    ctrlpp::mrac_controller<double> ctrl(cfg);

    ctrl.evaluate(vec1(0.0), vec1(1.0));
    auto xm1 = ctrl.x_model()[0];
    REQUIRE_THAT(xm1, WithinAbs(0.1, 1e-12));

    ctrl.evaluate(vec1(0.0), vec1(1.0));
    auto xm2 = ctrl.x_model()[0];
    REQUIRE_THAT(xm2, WithinAbs(0.9 * 0.1 + 0.1, 1e-12));
    REQUIRE(xm2 > xm1);
}
