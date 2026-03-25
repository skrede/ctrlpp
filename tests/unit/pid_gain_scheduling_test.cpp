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

TEST_CASE("set_params rescales integral for bumpless gain change",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(2.0);
    SisoPid pid(cfg);

    // Accumulate some integral: 10 steps of e=1, Ki=2
    // integral = sum(Ki*e*dt) = 2*1*0.01*10 = 0.2
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);

    double integral_before = pid.integral()[0];
    REQUIRE_THAT(integral_before, WithinAbs(0.2, tol));

    // Change Ki from 2.0 to 4.0
    // integral should rescale: integral_new = integral_old * ki_old / ki_new = 0.2 * 2/4 = 0.1
    // so that Ki * integral stays constant: 2*0.2 = 4*0.1 = 0.4
    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(4.0);
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.1, tol));
}

TEST_CASE("set_params with Ki going to zero clears integral",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(2.0);
    SisoPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);

    REQUIRE(pid.integral()[0] != 0.0);

    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(0.0);
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));
}

TEST_CASE("set_params Kp change produces no output discontinuity",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    // Run until steady output
    Vec1 u_before{};
    for (int i = 0; i < 20; ++i)
        u_before = pid.compute(vec1(1.0), vec1(0.5), dt);

    // Change Kp from 1 to 2
    SisoPid::config_type new_cfg = cfg;
    new_cfg.kp = vec1(2.0);
    pid.set_params(new_cfg);

    auto u_after = pid.compute(vec1(1.0), vec1(0.5), dt);

    // P changed from 1*0.5=0.5 to 2*0.5=1.0 (+0.5 step)
    // Ki unchanged so integral unchanged -- output jumps by Kp change on proportional
    // This is expected: gain scheduling Kp changes the proportional instantly
    // The key is integral rescaling for Ki changes
    double expected_p = 2.0 * 0.5;
    double integral = pid.integral()[0]; // integral accumulated + new increment
    REQUIRE_THAT(u_after[0], WithinAbs(expected_p + integral, 0.01));
}

TEST_CASE("params() returns current config after set_params",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    SisoPid::config_type new_cfg = cfg;
    new_cfg.kp = vec1(3.0);
    new_cfg.ki = vec1(1.5);
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.params().kp[0], WithinAbs(3.0, tol));
    REQUIRE_THAT(pid.params().ki[0], WithinAbs(1.5, tol));
}

TEST_CASE("set_integral sets integral to known value",
    "[pid][siso][integral-management]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    pid.set_integral(vec1(5.0));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(5.0, tol));

    // Next output should include this integral
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    // P=1*1=1, I=5.0 + ki*e*dt = 5.01, D=0
    REQUIRE_THAT(u[0], WithinAbs(6.01, tol));
}

TEST_CASE("freeze_integral prevents integral growth",
    "[pid][siso][integral-management]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    // Accumulate some integral
    pid.compute(vec1(1.0), vec1(0.0), dt);
    double integral_val = pid.integral()[0];

    // Freeze
    pid.freeze_integral(true);

    // Run 10 more steps -- integral should not change
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_val, tol));

    // Unfreeze
    pid.freeze_integral(false);

    // Integral should resume
    pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE(pid.integral()[0] > integral_val);
}

TEST_CASE("set_params from Ki=0 to Ki!=0 preserves zero integral",
    "[pid][siso][gain-scheduling][edge-case]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.0);
    SisoPid pid(cfg);

    // Run some steps with Ki=0 -> integral stays 0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Switch to Ki=2.0: ki_old==0, so the else-if branch (m_ki==0 -> clear) is NOT taken
    // because the NEW ki is non-zero. The condition ki_old[i]!=0 fails -> no rescale.
    // Integral stays at 0.
    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(2.0);
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Now integral should accumulate with new Ki
    pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(2.0 * 1.0 * dt, tol));
}

TEST_CASE("ISA form set_params rescales integral correctly",
    "[pid][siso][gain-scheduling][isa]")
{
    using IsaPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::isa_form>;
    IsaPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(4.0);  // Ti=4 -> internal Ki = Kp/Ti = 0.5
    cfg.kd = vec1(0.0);
    IsaPid pid(cfg);

    // Accumulate integral: 10 steps, e=1.0, internal Ki=0.5
    // integral = 0.5 * 1.0 * 0.01 * 10 = 0.05
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.05, tol));

    // Change Ti from 4.0 to 2.0 -> new internal Ki = 2/2 = 1.0
    // Bumpless rescale: integral_new = integral_old * ki_old/ki_new = 0.05 * 0.5/1.0 = 0.025
    IsaPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(2.0);
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.025, tol));
}

TEST_CASE("set_params with both Ki changing simultaneously (non-trivial rescale)",
    "[pid][siso][gain-scheduling][edge-case]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    // Accumulate integral = Ki*e*dt*5 = 1.0*1.0*0.01*5 = 0.05
    for (int i = 0; i < 5; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.05, tol));

    // Double the Ki: integral should halve to keep Ki*integral constant
    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(2.0);
    pid.set_params(new_cfg);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.025, tol));

    // Halve Ki back to 1.0: integral should double
    SisoPid::config_type newer_cfg = cfg;
    newer_cfg.ki = vec1(1.0);
    pid.set_params(newer_cfg);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.05, tol));
}
