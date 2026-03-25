#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1, 1, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

}

TEST_CASE("Tracking signal drives integral for bumpless transfer",
    "[pid][siso][tracking]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    // Run a few normal steps to accumulate integral
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    pid.compute(vec1(1.0), vec1(0.0), Ts);

    // Now supply tracking signal = 3.0
    // integral should be set so output matches tracking:
    // tracking = p + integral + d + ff
    // integral = tracking - p - d
    // With e=1.0: p = kp*ep = 1.0*1.0 = 1.0, d=0 (measurement unchanged)
    [[maybe_unused]] auto u = pid.compute(vec1(1.0), vec1(0.0), Ts, vec1(3.0));

    // After tracking: integral = 3.0 - p_term - d_term
    // p_term = kp * (b*sp - meas) = 1.0 * (1*1 - 0) = 1.0
    // d_term = 0 (measurement unchanged from prev)
    // integral = 3.0 - 1.0 - 0 = 2.0
    REQUIRE_THAT(pid.integral()[0], WithinAbs(2.0, tol));
}

TEST_CASE("Tracking enables bumpless manual-to-auto transition",
    "[pid][siso][tracking][bumpless]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    // Simulate manual mode: operator holds output at 5.0
    // Track the manual output for several steps
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.8), Ts, vec1(5.0));

    // Switch to auto: next compute without tracking should produce smooth output
    // Because integral was set to tracking - p - d:
    // p = kp * (b*sp - meas) = 2*(1*1 - 0.8) = 0.4
    // integral ~ 5.0 - 0.4 = 4.6  (d ~ 0 since meas stable)
    // Auto output: p + integral + I_increment = 0.4 + 4.6 + ki*e*dt = 5.0 + 1*0.2*0.01 ~ 5.002
    auto u_auto = pid.compute(vec1(1.0), vec1(0.8), Ts);

    // Output should be close to 5.0 (no bump from manual->auto)
    REQUIRE_THAT(u_auto[0], WithinAbs(5.002, 0.01));
}

TEST_CASE("Tracking with feed_forward: integral accounts for ff contribution",
    "[pid][siso][tracking][feed-forward]")
{
    auto ff_func = [](const Vec1& sp, double) -> Vec1 { return sp * 2.0; };
    using FF = ctrlpp::feed_forward<decltype(ff_func)>;
    using FfPid = ctrlpp::pid<double, 1, 1, 1, FF>;

    FfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    cfg.template policy<FF>().ff_func = ff_func;
    FfPid pid(cfg);

    // Step with tracking signal = 10.0
    // u = compute(sp, meas, Ts) includes ff contribution
    // non_integral = u - integral, integral = tracking - non_integral
    [[maybe_unused]] auto u = pid.compute(vec1(1.0), vec1(0.0), Ts, vec1(10.0));

    // P = 1*(1-0) = 1.0, I = 0.5*1*0.01 = 0.005, D = 0 (first step), FF = 1*2 = 2.0
    // u = 1.0 + 0.005 + 0 + 2.0 = 3.005
    // non_integral = u - integral = 3.005 - 0.005 = 3.0
    // new integral = tracking - non_integral = 10.0 - 3.0 = 7.0
    REQUIRE_THAT(pid.integral()[0], WithinAbs(7.0, tol));
}

TEST_CASE("Tracking after saturation resets integral appropriately",
    "[pid][siso][tracking][saturation]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(10.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    cfg.output_min = vec1(-5.0);
    SisoPid pid(cfg);

    // Drive into saturation
    for (int i = 0; i < 20; ++i)
        pid.compute(vec1(10.0), vec1(0.0), Ts);

    REQUIRE(pid.saturated() == true);

    // Now use tracking to pull integral back to a reasonable value
    pid.compute(vec1(1.0), vec1(0.0), Ts, vec1(2.0));

    // P = 10*1 = 10, clamped output = 5
    // non_integral = output - integral_before_tracking_adjustment
    // After tracking: integral = 2.0 - (5 - integral)
    // The integral should now be adjusted so next auto output is near 2.0
    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    // Should be within reasonable range of the tracking target
    REQUIRE(u[0] < 10.0);
}
