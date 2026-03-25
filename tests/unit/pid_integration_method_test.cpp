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

TEST_CASE("PI backward Euler accumulation", "[pid][siso][backward-euler]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    // Step 1: e=1.0, P=1.0, I=0.5*1.0*0.01=0.005
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.005, tol));

    // Step 2: e=1.0, P=1.0, I=0.005+0.005=0.01
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.01, tol));

    // Step 3: e=0.5, P=0.5, I=0.01+0.5*0.5*0.01=0.0125
    auto u3 = pid.compute(vec1(1.0), vec1(0.5), Ts);
    REQUIRE_THAT(u3[0], WithinAbs(0.5125, tol));
}

TEST_CASE("PI forward Euler uses previous error", "[pid][siso][forward-euler]")
{
    using FEPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::forward_euler>;
    FEPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    FEPid pid(cfg);

    // Step 1: prev_error = 0 (initialized), so I uses 0 * dt = 0
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    // P=1.0, I=0.5*0*0.01=0.0, total=1.0
    REQUIRE_THAT(u1[0], WithinAbs(1.0, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Step 2: prev_error = 1.0, I += 0.5*1.0*0.01 = 0.005
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(1.005, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.005, tol));
}

TEST_CASE("PI tustin uses trapezoidal average", "[pid][siso][tustin]")
{
    using TPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::tustin>;
    TPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    TPid pid(cfg);

    // Step 1: avg = (1.0 + 0.0) / 2 = 0.5, I = 0.5*0.5*0.01 = 0.0025
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(1.0025, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0025, tol));

    // Step 2: avg = (1.0 + 1.0) / 2 = 1.0, I += 0.5*1.0*0.01 = 0.005, total I = 0.0075
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(1.0075, tol));
}
