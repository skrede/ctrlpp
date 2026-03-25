#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

}

TEST_CASE("velocity_form P-only: delta_u = Kp*(e(k)-e(k-1))", "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    VPid pid(cfg);

    // Step 1: e=1.0, prev_e=0.0
    // delta_u = Kp*(1-0) + 0 + 0 = 2.0
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // Step 2: e=1.0, prev_e=1.0
    // delta_u = Kp*(1-1) = 0.0
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(0.0, tol));

    // Step 3: e=0.5, prev_e=1.0
    // delta_u = Kp*(0.5-1.0) = -1.0
    auto u3 = pid.compute(vec1(1.0), vec1(0.5), Ts);
    REQUIRE_THAT(u3[0], WithinAbs(-1.0, tol));
}

TEST_CASE("velocity_form PI: includes Ki*e*dt incremental term", "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    VPid pid(cfg);

    // Step 1: e=1.0, prev_e=0
    // dP = 1*(1-0) = 1.0, dI = 0.5*1*0.01 = 0.005
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    // Step 2: e=1.0, prev_e=1.0
    // dP = 1*(1-1) = 0, dI = 0.5*1*0.01 = 0.005
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(0.005, tol));

    // Step 3: e=0.5
    // dP = 1*(0.5-1) = -0.5, dI = 0.5*0.5*0.01 = 0.0025
    auto u3 = pid.compute(vec1(1.0), vec1(0.5), Ts);
    REQUIRE_THAT(u3[0], WithinAbs(-0.5 + 0.0025, tol));
}

TEST_CASE("velocity_form PID: full formula with second-order D difference",
    "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    cfg.kd = vec1(0.1);
    VPid pid(cfg);

    // Step 1: e(1)=1.0, e(0)=0, e(-1)=0
    // dP = 2*(1-0) = 2, dI = 0.5*1*0.01 = 0.005
    // dD = 0.1*(1 - 2*0 + 0)/0.01 = 0.1*100 = 10
    double du1 = 2.0 + 0.005 + 10.0;
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(du1, tol));

    // Step 2: e(2)=0.8, e(1)=1.0, e(0)=0
    // dP = 2*(0.8-1.0) = -0.4, dI = 0.5*0.8*0.01 = 0.004
    // dD = 0.1*(0.8 - 2*1.0 + 0)/0.01 = 0.1*(-1.2/0.01) = 0.1*(-120) = -12
    double du2 = -0.4 + 0.004 + (-12.0);
    auto u2 = pid.compute(vec1(1.0), vec1(0.2), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(du2, tol));

    // Step 3: e(3)=0.8, e(2)=0.8, e(1)=1.0
    // dP = 2*(0.8-0.8) = 0, dI = 0.5*0.8*0.01 = 0.004
    // dD = 0.1*(0.8 - 2*0.8 + 1.0)/0.01 = 0.1*(0.2/0.01) = 0.1*20 = 2
    double du3 = 0.0 + 0.004 + 2.0;
    auto u3 = pid.compute(vec1(1.0), vec1(0.2), Ts);
    REQUIRE_THAT(u3[0], WithinAbs(du3, tol));
}

TEST_CASE("velocity_form steady-state: delta_u converges to Ki*e*dt",
    "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(2.0);
    cfg.kd = vec1(0.1);
    VPid pid(cfg);

    // First few steps build up history, then with constant error:
    // dP = 0 (error unchanged), dD = 0 (e(k)-2e(k-1)+e(k-2)=0), dI = Ki*e*dt
    pid.compute(vec1(1.0), vec1(0.0), Ts);  // step 1
    pid.compute(vec1(1.0), vec1(0.0), Ts);  // step 2
    pid.compute(vec1(1.0), vec1(0.0), Ts);  // step 3

    // From step 3 onward, e is constant at 1.0
    // delta_u should be Ki*e*dt = 2*1*0.01 = 0.02
    auto u4 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u4[0], WithinAbs(0.02, tol));

    auto u5 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u5[0], WithinAbs(0.02, tol));
}

TEST_CASE("velocity_form + anti_windup compiles and runs (anti-windup is no-op)",
    "[pid][siso][velocity-form][anti-windup]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1,
        ctrlpp::velocity_form, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    VPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    VPid pid(cfg);

    // Should produce same result as without anti_windup
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    // dP = 1*(1-0) = 1, dI = 0.5*1*0.01 = 0.005
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));
}

TEST_CASE("velocity_form reset clears history", "[pid][siso][velocity-form][reset]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    VPid pid(cfg);

    // Build up history
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    pid.compute(vec1(1.0), vec1(0.0), Ts);

    pid.reset();

    // After reset: behaves as first step (prev_error=0, prev_prev_error=0)
    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
    // dP = 2*(1-0) = 2, dI = 0.5*1*0.01 = 0.005, dD = 0 (Kd=0)
    REQUIRE_THAT(u[0], WithinAbs(2.005, tol));
}

TEST_CASE("velocity_form with feed_forward adds ff delta to output",
    "[pid][siso][velocity-form][feed-forward]")
{
    auto ff_func = [](const Vec1& sp, double) -> Vec1 { return sp * 0.5; };
    using FF = ctrlpp::feed_forward<decltype(ff_func)>;
    using VFfPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form, FF>;

    VFfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.0);
    cfg.template policy<FF>().ff_func = ff_func;
    VFfPid pid(cfg);

    // Step 1: e=1, prev_e=0
    // dP = 1*(1-0) = 1.0, dI = 0, dD = 0
    // ff = sp * 0.5 = 0.5
    // delta_u = 1.0 + 0.5 = 1.5
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(1.5, tol));

    // Step 2: e=1, prev_e=1
    // dP = 1*(1-1) = 0, dI = 0, dD = 0
    // ff = sp * 0.5 = 0.5
    // delta_u = 0 + 0.5 = 0.5
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(0.5, tol));
}

TEST_CASE("velocity_form tracking signal is a no-op (integral not modified)",
    "[pid][siso][velocity-form][tracking]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    VPid pid(cfg);

    // Run a step with tracking; velocity form should ignore tracking signal
    // (the if constexpr branch for !velocity_form is not entered)
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts, vec1(99.0));
    // dP = 1*(1-0) = 1, dI = 0.5*1*0.01 = 0.005
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    // Integral should be zero (velocity form doesn't maintain position-form integral)
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));
}

TEST_CASE("velocity_form output clamping respects output_min and output_max",
    "[pid][siso][velocity-form][clamping]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(100.0);
    cfg.output_min = vec1(-0.5);
    cfg.output_max = vec1(0.5);
    VPid pid(cfg);

    // Step 1: dP = 100*(1-0) = 100 -> clamped to 0.5
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(0.5, tol));

    // Step 2: dP = 100*((-1)-1) = -200 -> clamped to -0.5
    auto u2 = pid.compute(vec1(0.0), vec1(1.0), Ts);
    REQUIRE_THAT(u2[0], WithinAbs(-0.5, tol));
}

TEST_CASE("velocity_form zero dt returns previous output",
    "[pid][siso][velocity-form][edge-case]")
{
    using VPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    VPid pid(cfg);

    auto u1 = pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // Zero dt returns previous output
    auto u2 = pid.compute(vec1(5.0), vec1(0.0), 0.0);
    REQUIRE_THAT(u2[0], WithinAbs(2.0, tol));
}
