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

TEST_CASE("Tracking signal drives integral for bumpless transfer",
    "[pid][siso][tracking]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    // Run a few normal steps to accumulate integral
    pid.compute(vec1(1.0), vec1(0.0), dt);
    pid.compute(vec1(1.0), vec1(0.0), dt);

    // Now supply tracking signal = 3.0
    // integral should be set so output matches tracking:
    // tracking = p + integral + d + ff
    // integral = tracking - p - d
    // With e=1.0: p = kp*ep = 1.0*1.0 = 1.0, d=0 (measurement unchanged)
    [[maybe_unused]] auto u = pid.compute(vec1(1.0), vec1(0.0), dt, vec1(3.0));

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
        pid.compute(vec1(1.0), vec1(0.8), dt, vec1(5.0));

    // Switch to auto: next compute without tracking should produce smooth output
    // Because integral was set to tracking - p - d:
    // p = kp * (b*sp - meas) = 2*(1*1 - 0.8) = 0.4
    // integral ~ 5.0 - 0.4 = 4.6  (d ~ 0 since meas stable)
    // Auto output: p + integral + I_increment = 0.4 + 4.6 + ki*e*dt = 5.0 + 1*0.2*0.01 ~ 5.002
    auto u_auto = pid.compute(vec1(1.0), vec1(0.8), dt);

    // Output should be close to 5.0 (no bump from manual->auto)
    REQUIRE_THAT(u_auto[0], WithinAbs(5.002, 0.01));
}
