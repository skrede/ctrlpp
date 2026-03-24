#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1, 1, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double dt = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

}

TEST_CASE("Output clamping", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(10.0);
    cfg.output_min = vec1(-5.0);
    cfg.output_max = vec1(5.0);
    SisoPid pid(cfg);

    // P = 10*1 = 10, clamped to 5
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(5.0, tol));

    // P = 10*(-1) = -10, clamped to -5
    auto u2 = pid.compute(vec1(0.0), vec1(1.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(-5.0, tol));
}

TEST_CASE("saturated() returns true when output is clamped", "[pid][siso][saturated]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(10.0);
    cfg.output_min = vec1(-5.0);
    cfg.output_max = vec1(5.0);
    SisoPid pid(cfg);

    // P = 10*1 = 10, clamped to 5 -> saturated
    pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE(pid.saturated() == true);

    // P = 10*0.1 = 1.0, not clamped -> not saturated
    pid.compute(vec1(0.1), vec1(0.0), dt);
    REQUIRE(pid.saturated() == false);
}

TEST_CASE("Bare pid with zero policies compiles", "[pid][siso][compile]")
{
    // This test verifies the zero-overhead case compiles
    using BarePid = ctrlpp::pid<double, 1, 1, 1>;
    BarePid::config_type cfg{};
    cfg.kp = vec1(1.0);
    BarePid pid(cfg);
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(1.0, tol));
}

TEST_CASE("Setpoint weighting b=0.5 modifies proportional term",
    "[pid][siso][setpoint-weighting]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.b = vec1(0.5);
    SisoPid pid(cfg);

    // ep = b*sp - meas = 0.5*1.0 - 0.0 = 0.5
    // P = 2.0 * 0.5 = 1.0
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(1.0, tol));
}

TEST_CASE("Setpoint weighting b=0 gives measurement-only proportional",
    "[pid][siso][setpoint-weighting]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.b = vec1(0.0);
    SisoPid pid(cfg);

    // ep = 0*sp - meas = -0.3
    // P = 2.0 * (-0.3) = -0.6
    auto u = pid.compute(vec1(1.0), vec1(0.3), dt);
    REQUIRE_THAT(u[0], WithinAbs(-0.6, tol));
}

TEST_CASE("Setpoint weighting does not affect integral (uses full error)",
    "[pid][siso][setpoint-weighting]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.b = vec1(0.5);
    SisoPid pid(cfg);

    pid.compute(vec1(2.0), vec1(0.0), dt);
    // Integral uses full error e = sp - meas = 2.0
    // I = Ki * e * dt = 1.0 * 2.0 * 0.01 = 0.02
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.02, tol));
}

TEST_CASE("Setpoint weighting c modifies derivative on error",
    "[pid][siso][setpoint-weighting][derivative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.kd = vec1(0.1);
    cfg.derivative_on_error = true;
    cfg.c = vec1(0.5);
    SisoPid pid(cfg);

    // Step 1: first step, D=0
    pid.compute(vec1(1.0), vec1(0.0), dt);

    // Step 2: sp jumps from 1.0 to 3.0, meas stays at 0.0
    // ed_curr = c*sp - meas = 0.5*3.0 - 0.0 = 1.5
    // ed_prev = c*prev_sp - prev_meas = 0.5*1.0 - 0.0 = 0.5
    // D = Kd * (ed_curr - ed_prev) / dt = 0.1 * (1.5 - 0.5) / 0.01 = 10.0
    // P = Kp * (b*sp - meas) = 1.0 * (1.0*3.0 - 0.0) = 3.0 (b default 1)
    // I = previous + Ki*e*dt (but ki=0 by default)
    auto u2 = pid.compute(vec1(3.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(3.0 + 10.0, tol));
}

TEST_CASE("Setpoint weighting c=0 with derivative on error uses only measurement change",
    "[pid][siso][setpoint-weighting][derivative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.kd = vec1(0.1);
    cfg.derivative_on_error = true;
    cfg.c = vec1(0.0);
    SisoPid pid(cfg);

    // Step 1
    pid.compute(vec1(1.0), vec1(0.0), dt);

    // Step 2: sp jumps from 1.0 to 3.0, meas stays 0
    // ed_curr = 0*3.0 - 0 = 0
    // ed_prev = 0*1.0 - 0 = 0
    // D = Kd * (0 - 0) / dt = 0
    // P = 1.0 * (1.0*3.0 - 0) = 3.0
    auto u2 = pid.compute(vec1(3.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(3.0, tol));
}

TEST_CASE("rate_limit constrains output change per step", "[pid][siso][rate-limit]")
{
    using RlPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::rate_limit>;
    RlPid::config_type cfg{};
    cfg.kp = vec1(100.0);
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {10.0};
    RlPid pid(cfg);

    // Step 1: prev_output = 0, raw = 100*1 = 100
    // max_delta = 10 * 0.01 = 0.1
    // rate limited: 0 + 0.1 = 0.1
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(0.1, tol));

    // Step 2: prev_output = 0.1, raw = 100*1 = 100
    // rate limited: 0.1 + 0.1 = 0.2
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(0.2, tol));
}

TEST_CASE("rate_limit applies before output clamp", "[pid][siso][rate-limit][pipeline-order]")
{
    using RlPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::rate_limit>;
    RlPid::config_type cfg{};
    cfg.kp = vec1(100.0);
    cfg.output_max = vec1(0.05);
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {10.0};
    RlPid pid(cfg);

    // Raw = 100, rate limited to 0.1, then clamped to 0.05
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(0.05, tol));
    REQUIRE(pid.saturated() == true);
}
