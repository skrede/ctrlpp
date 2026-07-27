#include "ctrlpp/control/l1.h"
#include "ctrlpp/types.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

auto vec1(double v) -> Vec1
{
    Vec1 r;
    r << v;
    return r;
}

// Plant: x[k+1] = 0.8*x + 0.5*u, DC gain = 2.5
// Predictor uses the same B as the plant for control effectiveness matching.
// A_m = 0.9 defines the desired closed-loop bandwidth.
// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. The rejection
// cases below do not use this helper: they assert the specific enumerator.
template <typename Controller, typename... Args>
auto make_controller(Args&&... args) -> Controller
{
    auto created = Controller::create(std::forward<Args>(args)...);
    REQUIRE(created.has_value());
    return *std::move(created);
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

TEST_CASE("l1 SISO step tracking", "[l1]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_controller<ctrlpp::l1_controller<double>>(cfg, 15.0, 100.0);

    auto r = vec1(1.0);
    double x_plant = 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        auto x = vec1(x_plant);
        auto u = ctrl.evaluate(x, r);
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }

    REQUIRE(std::abs(x_plant - 1.0) < 0.15);
}

TEST_CASE("l1 state predictor converges", "[l1]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_controller<ctrlpp::l1_controller<double>>(cfg, 15.0, 100.0);

    auto r = vec1(1.0);
    double x_plant = 0.0;

    for(int k = 0; k < 500; ++k)
    {
        auto x = vec1(x_plant);
        auto u = ctrl.evaluate(x, r);
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }

    REQUIRE(std::abs(ctrl.x_hat()[0] - x_plant) < 1.0);
}

TEST_CASE("l1 projection bounds respected", "[l1]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 100.0;
    cfg.theta_min << -2.0;
    cfg.theta_max << 2.0;

    auto ctrl = make_controller<ctrlpp::l1_controller<double>>(cfg, 5.0, 100.0);

    auto r = vec1(1.0);
    double x_plant = 0.0;

    for(int k = 0; k < 100; ++k)
    {
        auto x = vec1(x_plant);
        auto u = ctrl.evaluate(x, r);
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        REQUIRE(ctrl.sigma_hat()[0] >= -2.0);
        REQUIRE(ctrl.sigma_hat()[0] <= 2.0);
    }
}

TEST_CASE("l1 diagnostic accessors return expected types", "[l1]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_controller<ctrlpp::l1_controller<double>>(cfg, 5.0, 100.0);

    ctrl.evaluate(vec1(0.0), vec1(1.0));

    REQUIRE(ctrl.x_hat().size() == 1);
    REQUIRE(ctrl.sigma_hat().size() == 1);
    REQUIRE(ctrl.tracking_error().size() == 1);
    REQUIRE(ctrl.theta().size() == 1);
}

TEST_CASE("l1 reset restores initial state", "[l1]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_controller<ctrlpp::l1_controller<double>>(cfg, 5.0, 100.0);

    double x_plant = 0.0;
    for(int k = 0; k < 50; ++k)
    {
        auto u = ctrl.evaluate(vec1(x_plant), vec1(1.0));
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }

    REQUIRE(ctrl.x_hat()[0] != 0.0);

    ctrl.reset();

    REQUIRE_THAT(ctrl.x_hat()[0], WithinAbs(cfg.x_hat_0[0], 1e-15));
    REQUIRE_THAT(ctrl.sigma_hat()[0], WithinAbs(cfg.sigma_hat_0[0], 1e-15));
}

TEST_CASE("l1 direct constructor with vector_cascaded_biquad", "[l1]")
{
    auto cfg = make_siso_config();
    auto filter = ctrlpp::make_vector_butterworth<4, 1>(15.0, 100.0);
    REQUIRE(filter.has_value());
    auto ctrl = make_controller<
        ctrlpp::l1_controller<double, 1, 1, ctrlpp::vector_cascaded_biquad<double, 1, 2>>>(
        cfg, *std::move(filter));

    auto r = vec1(1.0);
    double x_plant = 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        auto x = vec1(x_plant);
        auto u = ctrl.evaluate(x, r);
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }

    REQUIRE(std::abs(x_plant - 1.0) < 0.15);
}

TEST_CASE("l1 MIMO 2x2 tracking", "[l1]")
{
    using Vec2 = ctrlpp::Vector<double, 2>;
    using Mat2 = ctrlpp::Matrix<double, 2, 2>;

    ctrlpp::l1_config<double, 2, 2> cfg{};
    cfg.predictor_model.A = (Mat2() << 0.9, 0.0, 0.0, 0.85).finished();
    cfg.predictor_model.B = (Mat2() << 0.5, 0.0, 0.0, 0.4).finished();
    cfg.predictor_model.C = Mat2::Identity();
    cfg.predictor_model.D = Mat2::Zero();
    cfg.gamma = 10.0 * Mat2::Identity();
    cfg.theta_min = (Vec2() << -5.0, -5.0).finished();
    cfg.theta_max = (Vec2() << 5.0, 5.0).finished();

    auto ctrl = make_controller<ctrlpp::l1_controller<double, 2, 2>>(cfg, 15.0, 100.0);

    Mat2 A_plant = (Mat2() << 0.8, 0.0, 0.0, 0.75).finished();
    Mat2 B_plant = (Mat2() << 0.5, 0.0, 0.0, 0.4).finished();
    Vec2 r = (Vec2() << 1.0, 0.5).finished();
    Vec2 x_plant = Vec2::Zero();

    for(int k = 0; k < 1000; ++k)
    {
        auto u = ctrl.evaluate(x_plant, r);
        x_plant = A_plant * x_plant + B_plant * u;
    }

    REQUIRE(std::abs(x_plant[0] - 1.0) < 0.2);
    REQUIRE(std::abs(x_plant[1] - 0.5) < 0.2);
}
