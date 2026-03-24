#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1, 1, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

}

TEST_CASE("MIMO NY=2 independent channels", "[pid][mimo]")
{
    using MimoPid = ctrlpp::pid<double, 1, 2, 2>;
    using Vec2 = ctrlpp::Vector<double, 2>;

    auto vec2 = [](double a, double b) { Vec2 v; v << a, b; return v; };

    MimoPid::config_type cfg{};
    cfg.kp = vec2(1.0, 3.0);
    cfg.ki = vec2(0.1, 0.2);
    MimoPid pid(cfg);

    Vec2 sp = vec2(1.0, 2.0);
    Vec2 meas = vec2(0.0, 0.0);

    auto u = pid.compute(sp, meas, 0.01);
    // Channel 0: P=1*1=1.0, I=0.1*1*0.01=0.001
    // Channel 1: P=3*2=6.0, I=0.2*2*0.01=0.004
    REQUIRE_THAT(u[0], WithinAbs(1.001, tol));
    REQUIRE_THAT(u[1], WithinAbs(6.004, tol));
}

TEST_CASE("Variable dt per step", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    // Step 1: dt=0.01, e=1, I=1*1*0.01=0.01
    [[maybe_unused]] auto u1 = pid.compute(vec1(1.0), vec1(0.0), 0.01);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.01, tol));

    // Step 2: dt=0.05, e=1, I=0.01+1*1*0.05=0.06
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), 0.05);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.06, tol));
    REQUIRE_THAT(u2[0], WithinAbs(1.06, tol));
}

TEST_CASE("MIMO performance metrics computed per channel",
    "[pid][mimo][perf-assessment]")
{
    using MimoPid = ctrlpp::pid<double, 1, 2, 2,
        ctrlpp::perf_assessment<ctrlpp::IAE>>;
    using Vec2 = ctrlpp::Vector<double, 2>;
    auto vec2 = [](double a, double b) { Vec2 v; v << a, b; return v; };
    MimoPid::config_type cfg{};
    cfg.kp = vec2(1.0, 1.0);
    MimoPid pid(cfg);

    // Channel 0: error=1.0, Channel 1: error=3.0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec2(1.0, 3.0), vec2(0.0, 0.0), 0.1);

    // IAE ch0 = 10*|1|*0.1 = 1.0, IAE ch1 = 10*|3|*0.1 = 3.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[1], WithinAbs(3.0, 1e-10));
}
