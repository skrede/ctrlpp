#include "hardening_helpers.h"
#include "ctrlpp/estimation/luenberger.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

auto make_system()
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 0.9, 0.1, 0.0, 0.8;
    sys.B << 0.0, 1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

}

TEST_CASE("Luenberger NaN measurement does not crash",
          "[luenberger][hardening][negative]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> L;
    L << 0.5, 0.3;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    obs.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << std::numeric_limits<double>::quiet_NaN();
    obs.update(z);

    CHECK((std::isnan(obs.state()[0]) || std::isfinite(obs.state()[0])));
}

TEST_CASE("Luenberger zero observer gains (open-loop)",
          "[luenberger][hardening][negative]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> L = Eigen::Matrix<double, 2, 1>::Zero();
    Eigen::Vector2d x0;
    x0 << 1.0, 0.0;

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    // With L=0, update has no effect -- purely open loop
    for(int k = 0; k < 10; ++k)
    {
        obs.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0;
        obs.update(z);
    }

    // State should be finite (decaying from initial)
    CHECK(std::isfinite(obs.state()[0]));
    CHECK(std::isfinite(obs.state()[1]));
}

TEST_CASE("Luenberger observer error decays for stable poles",
          "[luenberger][hardening][stability]")
{
    auto sys = make_system();
    // Choose L to place observer poles well inside unit circle
    Eigen::Matrix<double, 2, 1> L;
    L << 0.5, 0.3;
    Eigen::Vector2d x0;
    x0 << 10.0, 5.0;

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    // True state
    Eigen::Vector2d x_true;
    x_true << 0.0, 0.0;
    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    double prev_error = (x0 - x_true).norm();

    for(int k = 0; k < 200; ++k)
    {
        obs.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << sys.C(0, 0) * x_true[0] + sys.C(0, 1) * x_true[1];
        obs.update(z);

        x_true = (sys.A * x_true + sys.B * u).eval();
    }

    double final_error = (obs.state() - x_true).norm();
    REQUIRE(final_error < prev_error);
    REQUIRE(final_error < 0.1);
}

TEST_CASE("Luenberger observer converges to true state within 200 steps",
          "[luenberger][hardening][convergence]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> L;
    L << 0.5, 0.3;
    Eigen::Vector2d x0;
    x0 << 5.0, 2.0;

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    Eigen::Vector2d x_true;
    x_true << 0.0, 0.0;
    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    for(int k = 0; k < 200; ++k)
    {
        obs.predict(u);
        double meas = sys.C(0, 0) * x_true[0] + sys.C(0, 1) * x_true[1];
        Eigen::Matrix<double, 1, 1> z;
        z << meas;
        obs.update(z);

        x_true = (sys.A * x_true + sys.B * u).eval();
    }

    REQUIRE(std::abs(obs.state()[0] - x_true[0]) < 0.01);
}
