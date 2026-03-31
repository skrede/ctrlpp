#include "ctrlpp/control/mrac.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;

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

TEST_CASE("Dead-zone prevents parameter drift at steady state", "[mrac][dead_zone]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::dead_zone>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    cfg.robustification.threshold = 0.01;

    ctrl_t ctrl(cfg);

    double x = 0.0;
    double theta_x_5000 = 0.0;
    double theta_r_5000 = 0.0;

    for(int k = 0; k < 10000; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u[0];

        if(k == 4999)
        {
            theta_x_5000 = ctrl.theta_x()(0, 0);
            theta_r_5000 = ctrl.theta_r()(0, 0);
        }
    }

    REQUIRE(std::isfinite(ctrl.theta_x()(0, 0)));
    REQUIRE(std::isfinite(ctrl.theta_r()(0, 0)));
    CHECK(std::abs(ctrl.theta_x()(0, 0) - theta_x_5000) < 1e-6);
    CHECK(std::abs(ctrl.theta_r()(0, 0) - theta_r_5000) < 1e-6);
}

TEST_CASE("Dead-zone freezes adaptation below threshold", "[mrac][dead_zone]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::dead_zone>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    cfg.robustification.threshold = 100.0;

    ctrl_t ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 10; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u[0];
    }

    REQUIRE_THAT(ctrl.theta_x()(0, 0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ctrl.theta_r()(0, 0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("Sigma-modification adds leakage to parameters", "[mrac][sigma_modification]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::sigma_modification>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    cfg.theta_x_0 << 1.0;
    cfg.theta_r_0 << 1.0;
    cfg.robustification.sigma = 0.1;

    ctrl_t ctrl(cfg);

    ctrl.evaluate(vec1(0.0), vec1(0.0));

    REQUIRE_THAT(ctrl.theta_x()(0, 0), WithinAbs(0.9, 1e-10));
    REQUIRE_THAT(ctrl.theta_r()(0, 0), WithinAbs(0.9, 1e-10));
}

TEST_CASE("Sigma-modification parameters remain bounded over long horizon", "[mrac][sigma_modification]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::sigma_modification>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    cfg.robustification.sigma = 0.01;

    ctrl_t ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 10000; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u[0];
    }

    REQUIRE(std::isfinite(ctrl.theta_x()(0, 0)));
    REQUIRE(std::isfinite(ctrl.theta_r()(0, 0)));
    CHECK(std::abs(ctrl.theta_x()(0, 0)) < 1e6);
    CHECK(std::abs(ctrl.theta_r()(0, 0)) < 1e6);
}

TEST_CASE("e-modification scales leakage by error magnitude", "[mrac][e_modification]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::e_modification>;

    auto make_cfg = [](double x_model_init) {
        ctrl_t::config_type cfg{};
        cfg.reference_model = make_ref_model();
        cfg.gamma_x << 0.5;
        cfg.gamma_r << 0.5;
        cfg.theta_x_0 << 1.0;
        cfg.theta_r_0 << 1.0;
        cfg.robustification.delta = 0.1;
        cfg.x_model_0 << x_model_init;
        return cfg;
    };

    // Scenario 1: x_model_0 = 0 so e = x - x_model = 0 - x_model_after_propagation
    // With x=0 and r=0, x_model stays 0, so e=0 => no leakage
    ctrl_t ctrl_zero(make_cfg(0.0));
    ctrl_zero.evaluate(vec1(0.0), vec1(0.0));
    auto theta_x_zero_error = ctrl_zero.theta_x()(0, 0);

    // Scenario 2: x_model_0 = 10.0 so after propagation x_model = 0.9*10 = 9.0
    // e = 0 - 9.0 = -9.0, large error => significant leakage
    ctrl_t ctrl_large(make_cfg(10.0));
    ctrl_large.evaluate(vec1(0.0), vec1(0.0));
    auto theta_x_large_error = ctrl_large.theta_x()(0, 0);

    // With zero error, theta_x should remain at initial value (no adaptation, no leakage)
    CHECK_THAT(theta_x_zero_error, WithinAbs(1.0, 1e-10));

    // With large error, theta_x should have changed significantly due to leakage
    CHECK(std::abs(theta_x_large_error - 1.0) > std::abs(theta_x_zero_error - 1.0));
}

TEST_CASE("e-modification parameters remain bounded over long horizon", "[mrac][e_modification]")
{
    using ctrl_t = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::e_modification>;
    ctrl_t::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    cfg.robustification.delta = 0.01;

    ctrl_t ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 10000; ++k)
    {
        auto u = ctrl.evaluate(vec1(x), vec1(1.0));
        x = 0.8 * x + 0.5 * u[0];
    }

    REQUIRE(std::isfinite(ctrl.theta_x()(0, 0)));
    REQUIRE(std::isfinite(ctrl.theta_r()(0, 0)));
    CHECK(std::abs(ctrl.theta_x()(0, 0)) < 1e6);
    CHECK(std::abs(ctrl.theta_r()(0, 0)) < 1e6);
}
