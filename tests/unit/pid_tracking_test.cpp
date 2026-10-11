#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

template <typename Scalar, typename... Policies>
auto saturated_tracking(Scalar tracking) -> ctrlpp::pid<Scalar, 1, Policies...>
{
    using Pid = ctrlpp::pid<Scalar, 1, Policies...>;
    using Vec = ctrlpp::Vector<Scalar, 1>;
    typename Pid::config_type cfg{};
    cfg.kp = Vec::Constant(Scalar{1});
    cfg.ki = Vec::Constant(Scalar{1} / Scalar{10});
    cfg.output_min = Vec::Constant(Scalar{-10});
    cfg.output_max = Vec::Constant(Scalar{10});
    if constexpr(sizeof...(Policies) > 0)
        cfg.template policy<Policies...>().kb = {Scalar{1}};
    auto pid = ctrlpp::test::constructed(Pid::create(cfg));
    const Scalar dt = Scalar{1} / Scalar{100};
    REQUIRE(pid.compute(Vec::Constant(Scalar{20}), Vec::Constant(Scalar{0}), dt, Vec::Constant(tracking)).has_value());
    REQUIRE(pid.saturated());
    return pid;
}

}

TEST_CASE("Tracking signal drives integral for bumpless transfer",
    "[pid][siso][tracking]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    auto pid = ctrlpp::test::constructed(SisoPid::create(cfg));

    // Run a few normal steps to accumulate integral
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    // Now supply tracking signal = 3.0
    // integral should be set so output matches tracking:
    // tracking = p + integral + d + ff
    // integral = tracking - p - d
    // With e=1.0: p = kp*ep = 1.0*1.0 = 1.0, d=0 (measurement unchanged)
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts, vec1(3.0)).has_value());

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
    auto pid = ctrlpp::test::constructed(SisoPid::create(cfg));

    // Simulate manual mode: operator holds output at 5.0
    // Track the manual output for several steps
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(1.0), vec1(0.8), Ts, vec1(5.0)).has_value());

    // Switch to auto: next compute without tracking should produce smooth output
    // Because integral was set to tracking - p - d:
    // p = kp * (b*sp - meas) = 2*(1*1 - 0.8) = 0.4
    // integral ~ 5.0 - 0.4 = 4.6  (d ~ 0 since meas stable)
    // Auto output: p + integral + I_increment = 0.4 + 4.6 + ki*e*dt = 5.0 + 1*0.2*0.01 ~ 5.002
    auto u_auto = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.8), Ts));

    // Output should be close to 5.0 (no bump from manual->auto)
    REQUIRE_THAT(u_auto[0], WithinAbs(5.002, 0.01));
}

TEST_CASE("Tracking with feed_forward: integral accounts for ff contribution",
    "[pid][siso][tracking][feed-forward]")
{
    auto ff_func = [](const Vec1& sp, double) -> Vec1 { return sp * 2.0; };
    using FF = ctrlpp::feed_forward<decltype(ff_func)>;
    using FfPid = ctrlpp::pid<double, 1, FF>;

    FfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    cfg.template policy<FF>().ff_func = ff_func;
    auto pid = ctrlpp::test::constructed(FfPid::create(cfg));

    // Step with tracking signal = 10.0
    // integral = tracking - (P + D + FF)
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts, vec1(10.0)).has_value());

    // P = 1*(1-0) = 1.0, I = 0.5*1*0.01 = 0.005, D = 0 (first step), FF = 1*2 = 2.0
    // new integral = tracking - (P + D + FF) = 10.0 - 3.0 = 7.0
    REQUIRE_THAT(pid.integral()[0], WithinAbs(7.0, tol));
}

TEMPLATE_TEST_CASE("Tracking after saturation assigns the integrator from the unconstrained terms",
    "[pid][siso][tracking][saturation]", float, double)
{
    using Scalar = TestType;
    using Vec = ctrlpp::Vector<Scalar, 1>;
    const Scalar ki = Scalar{1} / Scalar{10};
    const Scalar dt = Scalar{1} / Scalar{100};
    const Scalar e = 20;
    auto pid = saturated_tracking<Scalar>(Scalar{5});
    REQUIRE(pid.integral()[0] == Scalar{-15});

    // The next automatic output is 5 before this cycle's own integral increment.
    const Scalar increment = (ki * e) * dt;
    const Scalar integral = Scalar{-15} + increment;
    const Scalar expected = (Scalar{20} + integral) + Scalar{0};
    auto u = ctrlpp::test::commanded(pid.compute(Vec::Constant(e), Vec::Constant(Scalar{0}), dt));
    REQUIRE(u[0] == expected);
}

TEMPLATE_TEST_CASE("Tracking clamps a target outside the output limits",
    "[pid][siso][tracking][saturation]", float, double)
{
    using Scalar = TestType;
    REQUIRE(saturated_tracking<Scalar>(Scalar{50}).integral()[0] == Scalar{-10});
    REQUIRE(saturated_tracking<Scalar>(Scalar{-50}).integral()[0] == Scalar{-30});
}

TEMPLATE_TEST_CASE("Tracking overrides the back-calculation step of the tracked cycle",
    "[pid][siso][tracking][saturation][backcalc]", float, double)
{
    using Scalar = TestType;
    auto pid = saturated_tracking<Scalar, ctrlpp::anti_windup<ctrlpp::back_calc>>(Scalar{5});
    REQUIRE(pid.integral()[0] == Scalar{-15});
}
