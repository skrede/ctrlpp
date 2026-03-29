#include "hardening_helpers.h"
#include "ctrlpp/control/l1.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

auto vec1(double v) -> Vec1
{
    Vec1 r;
    r << v;
    return r;
}

auto make_siso_config() -> ctrlpp::l1_config<double, 1, 1>
{
    ctrlpp::l1_config<double, 1, 1> cfg{};
    cfg.predictor_model.A << 0.9;
    cfg.predictor_model.B << 0.5;
    cfg.predictor_model.C << 1.0;
    cfg.predictor_model.D << 0.0;
    cfg.gamma << 10.0;
    cfg.theta_min << -10.0;
    cfg.theta_max << 10.0;
    return cfg;
}

}

TEST_CASE("L1 NaN state produces no crash", "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    ctrlpp::l1_controller<double> ctrl(cfg, 15.0, 100.0);

    auto u = ctrl.evaluate(ctrlpp::test::nan_vector<double, 1>(), vec1(1.0));
    CHECK((std::isnan(u[0]) || std::isfinite(u[0])));
}

TEST_CASE("L1 NaN reference produces no crash", "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    ctrlpp::l1_controller<double> ctrl(cfg, 15.0, 100.0);

    auto u = ctrl.evaluate(vec1(0.0), ctrlpp::test::nan_vector<double, 1>());
    CHECK((std::isnan(u[0]) || std::isfinite(u[0])));
}

TEST_CASE("L1 zero filter bandwidth produces finite output", "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    // Very low bandwidth filter -- should still produce finite output
    ctrlpp::l1_controller<double> ctrl(cfg, 0.1, 100.0);

    auto u = ctrl.evaluate(vec1(0.0), vec1(1.0));
    CHECK(std::isfinite(u[0]));
}

TEST_CASE("L1 known sigma_hat after one step", "[l1][hardening][precision]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1.0;
    ctrlpp::l1_controller<double> ctrl(cfg, 15.0, 100.0);

    // Step 1: x=0, r=1
    // x_hat = A_m*0 + B*(0+0) = 0
    // x_tilde = x_hat - x = 0 - 0 = 0
    // sigma_hat -= gamma * B^T * x_tilde = 0
    // sigma_hat clamped to [-10, 10] => 0
    auto u = ctrl.evaluate(vec1(0.0), vec1(1.0));
    REQUIRE_THAT(ctrl.sigma_hat()[0], WithinAbs(0.0, 1e-15));
    CHECK(std::isfinite(u[0]));
}

TEST_CASE("L1 tracks step reference with bounded transient", "[l1][hardening][convergence]")
{
    auto cfg = make_siso_config();
    ctrlpp::l1_controller<double> ctrl(cfg, 15.0, 100.0);

    double x_plant = 0.0;
    double max_overshoot = 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        auto u = ctrl.evaluate(vec1(x_plant), vec1(1.0));
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        if(x_plant > 1.0)
            max_overshoot = std::max(max_overshoot, x_plant - 1.0);
    }

    // Should converge close to reference
    REQUIRE(std::abs(x_plant - 1.0) < 0.15);
    // Overshoot should be bounded
    CHECK(max_overshoot < 1.0);
}

TEST_CASE("L1 high gamma low bandwidth bounded output", "[l1][hardening][robustness]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1e6;
    cfg.theta_min << -100.0;
    cfg.theta_max << 100.0;
    ctrlpp::l1_controller<double> ctrl(cfg, 2.0, 100.0);

    double x_plant = 0.0;
    bool all_finite = true;

    for(int k = 0; k < 500; ++k)
    {
        auto u = ctrl.evaluate(vec1(x_plant), vec1(1.0));
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        if(!std::isfinite(u[0]) || !std::isfinite(x_plant))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}

TEST_CASE("L1 projection clamps sigma_hat within bounds", "[l1][hardening][robustness]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1000.0;
    cfg.theta_min << -5.0;
    cfg.theta_max << 5.0;
    ctrlpp::l1_controller<double> ctrl(cfg, 15.0, 100.0);

    double x_plant = 0.0;

    for(int k = 0; k < 200; ++k)
    {
        auto u = ctrl.evaluate(vec1(x_plant), vec1(1.0));
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }

    CHECK(ctrl.sigma_hat()[0] >= -5.0 - 1e-10);
    CHECK(ctrl.sigma_hat()[0] <= 5.0 + 1e-10);
}
