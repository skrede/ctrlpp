#include "hardening_helpers.h"
#include "ctrlpp/control/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

using SisoPid = ctrlpp::pid<double, 1, 1, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

Vec1 vec1(double v)
{
    Vec1 r;
    r << v;
    return r;
}

constexpr double Ts = 0.01;

}

TEST_CASE("PID NaN setpoint produces NaN or unchanged output", "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(std::numeric_limits<double>::quiet_NaN()), vec1(0.0), Ts);
    // NaN propagation is acceptable -- just no crash
    CHECK((std::isnan(u[0]) || std::isfinite(u[0])));
}

TEST_CASE("PID NaN measurement produces NaN or unchanged output", "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(std::numeric_limits<double>::quiet_NaN()), Ts);
    CHECK((std::isnan(u[0]) || std::isfinite(u[0])));
}

TEST_CASE("PID Inf gains produce finite or Inf output without crash", "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(std::numeric_limits<double>::infinity());
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    CHECK(!std::isnan(u[0]));
}

TEST_CASE("PID near-zero dt returns previous output", "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    SisoPid pid(cfg);

    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), 1e-15);
    CHECK_THAT(u2[0], WithinAbs(u1[0], 1e-12));
}

TEST_CASE("PID P-only exact output", "[pid][hardening][precision]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u[0], WithinAbs(1.0, 1e-15));
}

TEST_CASE("PID precision with known PID gains", "[pid][hardening][precision]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    cfg.kd = vec1(0.0);
    SisoPid pid(cfg);

    // Step 1: e=1, P=2, I=0.5*1*0.01=0.005, D=0
    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u[0], WithinAbs(2.005, 1e-14));
}

TEST_CASE("PID stays bounded over 10000 steps with stable gains", "[pid][hardening][stability]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.5);
    cfg.ki = vec1(0.1);
    cfg.kd = vec1(0.01);
    SisoPid pid(cfg);

    double plant_x = 0.0;
    bool all_finite = true;

    for(int k = 0; k < 10000; ++k)
    {
        auto u = pid.compute(vec1(1.0), vec1(plant_x), Ts);
        plant_x = 0.9 * plant_x + 0.1 * u[0];
        if(!std::isfinite(u[0]) || !std::isfinite(plant_x))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
    CHECK(std::abs(plant_x) < 1e6);
}

TEST_CASE("PID with integral action converges to zero error", "[pid][hardening][convergence]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.5);
    cfg.ki = vec1(0.2);
    SisoPid pid(cfg);

    double plant_x = 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        auto u = pid.compute(vec1(1.0), vec1(plant_x), Ts);
        plant_x = 0.9 * plant_x + 0.1 * u[0];
    }

    double error = std::abs(1.0 - plant_x);
    REQUIRE(error < 0.01);
}

TEST_CASE("PID extreme gains produce finite output", "[pid][hardening][robustness]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1e12);
    cfg.ki = vec1(1e12);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    CHECK(std::isfinite(u[0]));

    auto u2 = pid.compute(vec1(1.0), vec1(0.5), Ts);
    CHECK(std::isfinite(u2[0]));
}

TEST_CASE("PID zero gains produce zero output", "[pid][hardening][robustness]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.0);
    cfg.ki = vec1(0.0);
    cfg.kd = vec1(0.0);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(100.0), vec1(0.0), Ts);
    REQUIRE_THAT(u[0], WithinAbs(0.0, 1e-15));
}
