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

TEST_CASE("Without anti_windup: integral winds up unboundedly", "[pid][siso][windup]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    SisoPid pid(cfg);

    // Constant error of 1.0 for many steps
    for (int i = 0; i < 1000; ++i)
        pid.compute(vec1(10.0), vec1(0.0), dt);

    // Integral should be large (no anti-windup to stop it)
    // I = Ki * e * dt * 1000 = 1 * 10 * 0.01 * 1000 ~ 100
    REQUIRE_THAT(pid.integral()[0], WithinAbs(100.0, 1e-6));
}

TEST_CASE("back_calc anti-windup limits integral growth during saturation",
    "[pid][siso][anti-windup][backcalc]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::back_calc>;
    using AwPid = ctrlpp::pid<double, 1, 1, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    cfg.template policy<AW>().kb = {1.0};
    AwPid pid(cfg);

    // Run many steps with large error to cause saturation
    for (int i = 0; i < 1000; ++i)
        pid.compute(vec1(10.0), vec1(0.0), dt);

    // Integral should be bounded (back-calc feedback limits growth)
    // Without anti-windup it would be 100.0
    REQUIRE(pid.integral()[0] < 50.0);
}

TEST_CASE("back_calc default Kb auto-computation", "[pid][siso][anti-windup][backcalc][auto-kb]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::back_calc>;

    SECTION("Kb = sqrt(Ki*Kd) when Kd != 0") {
        using AwPid = ctrlpp::pid<double, 1, 1, 1, AW>;
        AwPid::config_type cfg{};
        cfg.kp = vec1(2.0);
        cfg.ki = vec1(4.0);
        cfg.kd = vec1(1.0);
        cfg.output_max = vec1(5.0);
        // kb left at 0 -> auto-compute = sqrt(4*1) = 2
        AwPid pid(cfg);

        // Run until saturation
        for (int i = 0; i < 100; ++i)
            pid.compute(vec1(10.0), vec1(0.0), dt);

        // Verify anti-windup is active (integral bounded)
        REQUIRE(pid.integral()[0] < 50.0);
    }

    SECTION("Kb = Ki when Kd == 0") {
        using AwPid = ctrlpp::pid<double, 1, 1, 1, AW>;
        AwPid::config_type cfg{};
        cfg.kp = vec1(1.0);
        cfg.ki = vec1(2.0);
        // kd left at 0
        cfg.output_max = vec1(5.0);
        AwPid pid(cfg);

        for (int i = 0; i < 100; ++i)
            pid.compute(vec1(10.0), vec1(0.0), dt);

        REQUIRE(pid.integral()[0] < 50.0);
    }
}

TEST_CASE("clamping anti-windup freezes integral during saturation",
    "[pid][siso][anti-windup][clamping]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::clamping>;
    using AwPid = ctrlpp::pid<double, 1, 1, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    AwPid pid(cfg);

    // Step until saturation
    double integral_at_saturation = 0.0;
    bool found_saturation = false;
    for (int i = 0; i < 100; ++i) {
        pid.compute(vec1(10.0), vec1(0.0), dt);
        if (pid.saturated() && !found_saturation) {
            found_saturation = true;
            integral_at_saturation = pid.integral()[0];
        }
    }
    REQUIRE(found_saturation);

    // After many more steps during saturation, integral should be frozen
    // (error > 0 and integral > 0 during saturation -> undo increment)
    for (int i = 0; i < 100; ++i)
        pid.compute(vec1(10.0), vec1(0.0), dt);

    // Integral should have stayed near the saturation point
    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_at_saturation, tol));
}

TEST_CASE("conditional_integration freezes integral when error exceeds threshold",
    "[pid][siso][anti-windup][conditional]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::conditional_integration>;
    using AwPid = ctrlpp::pid<double, 1, 1, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.template policy<AW>().error_threshold = {2.0};
    AwPid pid(cfg);

    // Error = 5 > threshold 2 -> integral should not accumulate
    for (int i = 0; i < 100; ++i)
        pid.compute(vec1(5.0), vec1(0.0), dt);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Error = 1 < threshold 2 -> integral should accumulate
    pid.reset();
    pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(1.0 * 1.0 * dt, tol));
}
