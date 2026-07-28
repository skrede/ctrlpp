#include "hardening_helpers.h"
#include "ctrlpp/control/mrac.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;

namespace {

using Vec2 = ctrlpp::Vector<double, 2>;
using Mat2 = ctrlpp::Matrix<double, 2, 2>;

auto make_mimo_ref_model() -> ctrlpp::discrete_state_space<double, 2, 2, 2>
{
    ctrlpp::discrete_state_space<double, 2, 2, 2> ref{};
    ref.A << 0.9, 0.0,
             0.0, 0.85;
    ref.B << 0.1, 0.0,
             0.0, 0.15;
    ref.C = Mat2::Identity();
    ref.D = Mat2::Zero();
    return ref;
}

auto make_mimo_config() -> ctrlpp::mrac_config<double, 2, 2>
{
    ctrlpp::mrac_config<double, 2, 2> cfg{};
    cfg.reference_model = make_mimo_ref_model();
    cfg.gamma_x = 0.3 * Mat2::Identity();
    cfg.gamma_r = 0.3 * Mat2::Identity();
    return cfg;
}

}

TEST_CASE("MIMO MRAC construction and accessors", "[mrac][mimo]")
{
    auto cfg = make_mimo_config();
    ctrlpp::mrac_controller<double, 2, 2> ctrl(cfg);

    REQUIRE(ctrl.theta_x().norm() < 1e-15);
    REQUIRE(ctrl.theta_r().norm() < 1e-15);
    REQUIRE_THAT(ctrl.x_model()[0], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.x_model()[1], WithinAbs(0.0, 1e-15));
}

TEST_CASE("MIMO MRAC tracks 2D step reference", "[mrac][mimo]")
{
    auto cfg = make_mimo_config();
    ctrlpp::mrac_controller<double, 2, 2> ctrl(cfg);

    Mat2 A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    Mat2 B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    Vec2 x = Vec2::Zero();
    Vec2 r;
    r << 1.0, 0.5;

    for(int k = 0; k < 500; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(x, r));
        x = A_p * x + B_p * u;
    }

    auto x_m = ctrl.x_model();
    REQUIRE((x - x_m).norm() < 0.1 * x_m.norm());
}

TEST_CASE("MIMO MRAC evaluate returns 2D control vector", "[mrac][mimo]")
{
    auto cfg = make_mimo_config();
    ctrlpp::mrac_controller<double, 2, 2> ctrl(cfg);

    Vec2 x = Vec2::Zero();
    Vec2 r;
    r << 1.0, 0.5;

    auto u = ctrlpp::test::commanded(ctrl.evaluate(x, r));

    REQUIRE(std::isfinite(u[0]));
    REQUIRE(std::isfinite(u[1]));
}

TEST_CASE("MIMO dead-zone prevents adaptation below threshold", "[mrac][mimo]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 2, 2, ctrlpp::dead_zone>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_mimo_ref_model();
    cfg.gamma_x = 0.3 * Mat2::Identity();
    cfg.gamma_r = 0.3 * Mat2::Identity();
    cfg.robustification.threshold = 100.0;

    ctrl_t ctrl(cfg);

    Mat2 A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    Mat2 B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    Vec2 x = Vec2::Zero();
    Vec2 r;
    r << 1.0, 0.5;

    for(int k = 0; k < 10; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(x, r));
        x = A_p * x + B_p * u;
    }

    REQUIRE(ctrl.theta_x().norm() < 1e-14);
    REQUIRE(ctrl.theta_r().norm() < 1e-14);
}

TEST_CASE("MIMO sigma-modification bounds parameters", "[mrac][mimo]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 2, 2, ctrlpp::sigma_modification>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_mimo_ref_model();
    cfg.gamma_x = 0.3 * Mat2::Identity();
    cfg.gamma_r = 0.3 * Mat2::Identity();
    cfg.robustification.sigma = 0.01;

    ctrl_t ctrl(cfg);

    Mat2 A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    Mat2 B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    Vec2 x = Vec2::Zero();
    Vec2 r;
    r << 1.0, 0.5;

    for(int k = 0; k < 10000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(x, r));
        x = A_p * x + B_p * u;
    }

    REQUIRE(std::isfinite(ctrl.theta_x().norm()));
    REQUIRE(std::isfinite(ctrl.theta_r().norm()));
    CHECK(ctrl.theta_x().norm() < 1e6);
    CHECK(ctrl.theta_r().norm() < 1e6);
}

TEST_CASE("MIMO e-modification bounds parameters", "[mrac][mimo]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 2, 2, ctrlpp::e_modification>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_mimo_ref_model();
    cfg.gamma_x = 0.3 * Mat2::Identity();
    cfg.gamma_r = 0.3 * Mat2::Identity();
    cfg.robustification.delta = 0.01;

    ctrl_t ctrl(cfg);

    Mat2 A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    Mat2 B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    Vec2 x = Vec2::Zero();
    Vec2 r;
    r << 1.0, 0.5;

    for(int k = 0; k < 10000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(x, r));
        x = A_p * x + B_p * u;
    }

    REQUIRE(std::isfinite(ctrl.theta_x().norm()));
    REQUIRE(std::isfinite(ctrl.theta_r().norm()));
    CHECK(ctrl.theta_x().norm() < 1e6);
    CHECK(ctrl.theta_r().norm() < 1e6);
}

TEST_CASE("MIMO MRAC reset restores initial state", "[mrac][mimo]")
{
    auto cfg = make_mimo_config();
    ctrlpp::mrac_controller<double, 2, 2> ctrl(cfg);

    Mat2 A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    Mat2 B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    Vec2 x = Vec2::Zero();
    Vec2 r;
    r << 1.0, 0.5;

    for(int k = 0; k < 50; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(x, r));
        x = A_p * x + B_p * u;
    }

    REQUIRE(ctrl.theta_x().norm() > 1e-10);

    ctrl.reset();

    REQUIRE(ctrl.theta_x().norm() < 1e-15);
    REQUIRE(ctrl.theta_r().norm() < 1e-15);
    REQUIRE_THAT(ctrl.x_model()[0], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.x_model()[1], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.tracking_error()[0], WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.tracking_error()[1], WithinAbs(0.0, 1e-15));
}
