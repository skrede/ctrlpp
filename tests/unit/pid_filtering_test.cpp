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

TEST_CASE("setpoint_filter smooths step setpoint", "[pid][siso][setpoint-filter]")
{
    using SpfPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::setpoint_filter>;
    SpfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.template policy<ctrlpp::setpoint_filter>().tf = {0.1};
    SpfPid pid(cfg);

    // alpha = 0.1 / (0.1 + 0.01) = 10/11 ~ 0.9091
    double alpha = 0.1 / (0.1 + Ts);

    // Step 1: filtered_sp = alpha*0 + (1-alpha)*1.0 = (1-alpha)
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    double fsp1 = (1.0 - alpha) * 1.0;
    REQUIRE_THAT(u1[0], WithinAbs(fsp1, tol));  // P = 1.0 * (fsp1 - 0)

    // Step 10: filtered_sp approaches 1.0 gradually
    double fsp = fsp1;
    for (int i = 1; i < 10; ++i) {
        fsp = alpha * fsp + (1.0 - alpha) * 1.0;
        pid.compute(vec1(1.0), vec1(0.0), Ts);
    }
    // After 10 steps, filtered_sp should be closer to 1.0 but not yet there
    REQUIRE(fsp > 0.5);
    REQUIRE(fsp < 1.0);
}

TEST_CASE("pv_filter smooths step measurement", "[pid][siso][pv-filter]")
{
    using PvfPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::pv_filter>;
    PvfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.template policy<ctrlpp::pv_filter>().tf = {0.05};
    PvfPid pid(cfg);

    double alpha = 0.05 / (0.05 + Ts);

    // Step 1: filtered_meas = (1-alpha)*1.0
    auto u1 = pid.compute(vec1(0.0), vec1(1.0), Ts);
    double fm1 = (1.0 - alpha) * 1.0;
    // e = 0 - fm1 = -fm1, P = -fm1
    REQUIRE_THAT(u1[0], WithinAbs(-fm1, tol));

    // After several steps, filtered_meas approaches 1.0
    double fm = fm1;
    for (int i = 1; i < 20; ++i) {
        fm = alpha * fm + (1.0 - alpha) * 1.0;
        pid.compute(vec1(0.0), vec1(1.0), Ts);
    }
    REQUIRE(fm > 0.9);
    REQUIRE(fm < 1.0);
}

TEST_CASE("feed_forward adds callable output to control signal", "[pid][siso][feed-forward]")
{
    // Feed-forward: ff(sp, Ts) = sp (identity)
    auto ff_func = [](const Vec1& sp, double) -> Vec1 { return sp; };
    using FF = ctrlpp::feed_forward<decltype(ff_func)>;
    using FfPid = ctrlpp::pid<double, 1, 1, 1, FF>;

    FfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.template policy<FF>().ff_func = ff_func;
    FfPid pid(cfg);

    // P = 1.0 * (1.0 - 0.0) = 1.0, FF = sp = 1.0, total = 2.0
    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u[0], WithinAbs(2.0, tol));
}

TEST_CASE("setpoint_filter + pv_filter combined", "[pid][siso][filter-combo]")
{
    using ComboPid = ctrlpp::pid<double, 1, 1, 1,
        ctrlpp::setpoint_filter, ctrlpp::pv_filter>;
    ComboPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.template policy<ctrlpp::setpoint_filter>().tf = {0.1};
    cfg.template policy<ctrlpp::pv_filter>().tf = {0.05};
    ComboPid pid(cfg);

    double alpha_sp = 0.1 / (0.1 + Ts);
    double alpha_pv = 0.05 / (0.05 + Ts);

    // Step 1: sp=1.0, meas=0.5
    double fsp = (1.0 - alpha_sp) * 1.0;
    double fm = (1.0 - alpha_pv) * 0.5;
    auto u = pid.compute(vec1(1.0), vec1(0.5), Ts);
    // e = fsp - fm, P = 1.0 * e
    REQUIRE_THAT(u[0], WithinAbs(fsp - fm, tol));
}

TEST_CASE("No filter policies: output unchanged from baseline", "[pid][siso][no-filter]")
{
    // Bare PID (no filter policies) should match Plan 01 baseline exactly
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    // P=2*1=2, I=0.5*1*0.01=0.005
    REQUIRE_THAT(u[0], WithinAbs(2.005, tol));
}

TEST_CASE("setpoint_filter reset clears filter state", "[pid][siso][setpoint-filter][reset]")
{
    using SpfPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::setpoint_filter>;
    SpfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.template policy<ctrlpp::setpoint_filter>().tf = {0.1};
    SpfPid pid(cfg);

    double alpha = 0.1 / (0.1 + Ts);

    // Run a few steps to build up filter state
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    pid.compute(vec1(1.0), vec1(0.0), Ts);

    pid.reset();

    // After reset, filter state is cleared -- first step behaves as if fresh
    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    double fsp = (1.0 - alpha) * 1.0;
    REQUIRE_THAT(u[0], WithinAbs(fsp, tol));
}
