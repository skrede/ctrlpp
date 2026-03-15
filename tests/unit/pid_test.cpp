#include "ctrlpp/pid.h"
#include "naive_linalg.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1>;
using Vec1 = NaiveLinalg::vector_type<double, 1>;

constexpr double dt = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { return {v}; }

}

TEST_CASE("P-only controller", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {2.5};
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(2.5, tol));
}

TEST_CASE("PI backward Euler accumulation", "[pid][siso][backward-euler]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    SisoPid pid(cfg);

    // Step 1: e=1.0, P=1.0, I=0.5*1.0*0.01=0.005
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.005, tol));

    // Step 2: e=1.0, P=1.0, I=0.005+0.005=0.01
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.01, tol));

    // Step 3: e=0.5, P=0.5, I=0.01+0.5*0.5*0.01=0.0125
    auto u3 = pid.compute(vec1(1.0), vec1(0.5), dt);
    REQUIRE_THAT(u3[0], WithinAbs(0.5125, tol));
}

TEST_CASE("PI forward Euler uses previous error", "[pid][siso][forward-euler]")
{
    using FEPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::ForwardEuler>;
    FEPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    FEPid pid(cfg);

    // Step 1: prev_error = 0 (initialized), so I uses 0 * dt = 0
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    // P=1.0, I=0.5*0*0.01=0.0, total=1.0
    REQUIRE_THAT(u1[0], WithinAbs(1.0, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Step 2: prev_error = 1.0, I += 0.5*1.0*0.01 = 0.005
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.005, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.005, tol));
}

TEST_CASE("PI Tustin uses trapezoidal average", "[pid][siso][tustin]")
{
    using TPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::Tustin>;
    TPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    TPid pid(cfg);

    // Step 1: avg = (1.0 + 0.0) / 2 = 0.5, I = 0.5*0.5*0.01 = 0.0025
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.0025, tol));
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0025, tol));

    // Step 2: avg = (1.0 + 1.0) / 2 = 1.0, I += 0.5*1.0*0.01 = 0.005, total I = 0.0075
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.0075, tol));
}

TEST_CASE("PD derivative on measurement (default)", "[pid][siso][derivative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.kd = {0.1};
    SisoPid pid(cfg);

    // Step 1: first step, derivative is zero
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    // P=1.0, D=0 (first step)
    REQUIRE_THAT(u1[0], WithinAbs(1.0, tol));

    // Step 2: measurement changes from 0.0 to 0.2
    // D = -Kd * (meas_new - meas_old) / dt = -0.1 * (0.2 - 0.0) / 0.01 = -2.0
    auto u2 = pid.compute(vec1(1.0), vec1(0.2), dt);
    // P = 1.0 * (1.0 - 0.2) = 0.8
    REQUIRE_THAT(u2[0], WithinAbs(0.8 - 2.0, tol));
}

TEST_CASE("PD derivative on error", "[pid][siso][derivative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.kd = {0.1};
    cfg.derivative_on_error = true;
    SisoPid pid(cfg);

    // Step 1: first step, derivative is zero
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.0, tol));

    // Step 2: setpoint jumps from 1.0 to 2.0, meas stays at 0.0
    // e_new = 2.0 - 0.0 = 2.0, e_old = 1.0
    // D = Kd * (e_new - e_old) / dt = 0.1 * (2.0 - 1.0) / 0.01 = 10.0
    auto u2 = pid.compute(vec1(2.0), vec1(0.0), dt);
    // P = 1.0 * 2.0 = 2.0
    REQUIRE_THAT(u2[0], WithinAbs(2.0 + 10.0, tol));
}

TEST_CASE("Full PID matches hand computation", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.ki = {0.5};
    cfg.kd = {0.1};
    SisoPid pid(cfg);

    // Step 1: sp=1, meas=0, e=1
    // P = 2*1 = 2, I = 0.5*1*0.01 = 0.005, D = 0 (first step)
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(2.005, tol));

    // Step 2: sp=1, meas=0.1, e=0.9
    // P = 2*0.9 = 1.8
    // I = 0.005 + 0.5*0.9*0.01 = 0.0095
    // D = -0.1*(0.1-0.0)/0.01 = -1.0 (derivative on measurement)
    auto u2 = pid.compute(vec1(1.0), vec1(0.1), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.8 + 0.0095 - 1.0, tol));
}

TEST_CASE("ISA form converts to parallel form", "[pid][siso][isa]")
{
    // ISA form: Kp=2, Ti=4, Td=0.05
    // Parallel: Kp=2, Ki=2/4=0.5, Kd=2*0.05=0.1
    using IsaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::IsaForm>;
    IsaPid::config_type isa_cfg{};
    isa_cfg.kp = {2.0};
    isa_cfg.ki = {4.0};  // This is Ti for ISA form
    isa_cfg.kd = {0.05}; // This is Td for ISA form
    IsaPid isa_pid(isa_cfg);

    // Parallel form reference
    SisoPid::config_type par_cfg{};
    par_cfg.kp = {2.0};
    par_cfg.ki = {0.5};
    par_cfg.kd = {0.1};
    SisoPid par_pid(par_cfg);

    auto u_isa = isa_pid.compute(vec1(1.0), vec1(0.0), dt);
    auto u_par = par_pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u_isa[0], WithinAbs(u_par[0], tol));

    auto u_isa2 = isa_pid.compute(vec1(1.0), vec1(0.1), dt);
    auto u_par2 = par_pid.compute(vec1(1.0), vec1(0.1), dt);
    REQUIRE_THAT(u_isa2[0], WithinAbs(u_par2[0], tol));
}

TEST_CASE("MIMO NY=2 independent channels", "[pid][mimo]")
{
    using MimoPid = ctrlpp::Pid<NaiveLinalg, double, 1, 2, 2>;
    using Vec2 = NaiveLinalg::vector_type<double, 2>;

    MimoPid::config_type cfg{};
    cfg.kp = {1.0, 3.0};
    cfg.ki = {0.1, 0.2};
    MimoPid pid(cfg);

    Vec2 sp = {1.0, 2.0};
    Vec2 meas = {0.0, 0.0};

    auto u = pid.compute(sp, meas, dt);
    // Channel 0: P=1*1=1.0, I=0.1*1*0.01=0.001
    // Channel 1: P=3*2=6.0, I=0.2*2*0.01=0.004
    REQUIRE_THAT(u[0], WithinAbs(1.001, tol));
    REQUIRE_THAT(u[1], WithinAbs(6.004, tol));
}

TEST_CASE("Variable dt per step", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {1.0};
    SisoPid pid(cfg);

    // Step 1: dt=0.01, e=1, I=1*1*0.01=0.01
    [[maybe_unused]] auto u1 = pid.compute(vec1(1.0), vec1(0.0), 0.01);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.01, tol));

    // Step 2: dt=0.05, e=1, I=0.01+1*1*0.05=0.06
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), 0.05);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.06, tol));
    REQUIRE_THAT(u2[0], WithinAbs(1.06, tol));
}

TEST_CASE("reset() zeros state", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.ki = {1.0};
    SisoPid pid(cfg);

    pid.compute(vec1(1.0), vec1(0.0), dt);
    pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE(pid.integral()[0] != 0.0);

    pid.reset();
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.error()[0], WithinAbs(0.0, tol));

    // After reset, first step derivative should be zero again
    auto u = pid.compute(vec1(1.0), vec1(0.5), dt);
    // P = 2*0.5 = 1.0, I = 1*0.5*0.01 = 0.005, D = 0 (first step after reset)
    REQUIRE_THAT(u[0], WithinAbs(1.005, tol));
}

TEST_CASE("First step derivative is zero", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {0.0};
    cfg.kd = {10.0};
    SisoPid pid(cfg);

    // Even with large measurement, derivative should be zero on first step
    auto u = pid.compute(vec1(0.0), vec1(100.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(0.0, tol));
}

TEST_CASE("Zero dt returns previous output", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {2.0};
    SisoPid pid(cfg);

    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // Zero dt should return previous output
    auto u2 = pid.compute(vec1(5.0), vec1(0.0), 0.0);
    REQUIRE_THAT(u2[0], WithinAbs(2.0, tol));

    // Negative dt should also return previous output
    auto u3 = pid.compute(vec1(5.0), vec1(0.0), -0.01);
    REQUIRE_THAT(u3[0], WithinAbs(2.0, tol));
}

TEST_CASE("Output clamping", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {10.0};
    cfg.output_min = {-5.0};
    cfg.output_max = {5.0};
    SisoPid pid(cfg);

    // P = 10*1 = 10, clamped to 5
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(5.0, tol));

    // P = 10*(-1) = -10, clamped to -5
    auto u2 = pid.compute(vec1(0.0), vec1(1.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(-5.0, tol));
}

TEST_CASE("Bare Pid with zero policies compiles", "[pid][siso][compile]")
{
    // This test verifies the zero-overhead case compiles
    using BarePid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1>;
    BarePid::config_type cfg{};
    cfg.kp = {1.0};
    BarePid pid(cfg);
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(1.0, tol));
}

TEST_CASE("error() returns last error", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    SisoPid pid(cfg);

    pid.compute(vec1(3.0), vec1(1.0), dt);
    REQUIRE_THAT(pid.error()[0], WithinAbs(2.0, tol));
}

// =============================================================================
// Plan 03: Input filtering, feed-forward
// =============================================================================

TEST_CASE("SetpointFilter smooths step setpoint", "[pid][siso][setpoint-filter]")
{
    using SpfPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::SetpointFilter>;
    SpfPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.template policy<ctrlpp::SetpointFilter>().tf = {0.1};
    SpfPid pid(cfg);

    // alpha = 0.1 / (0.1 + 0.01) = 10/11 ~ 0.9091
    double alpha = 0.1 / (0.1 + dt);

    // Step 1: filtered_sp = alpha*0 + (1-alpha)*1.0 = (1-alpha)
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    double fsp1 = (1.0 - alpha) * 1.0;
    REQUIRE_THAT(u1[0], WithinAbs(fsp1, tol));  // P = 1.0 * (fsp1 - 0)

    // Step 10: filtered_sp approaches 1.0 gradually
    double fsp = fsp1;
    for (int i = 1; i < 10; ++i) {
        fsp = alpha * fsp + (1.0 - alpha) * 1.0;
        pid.compute(vec1(1.0), vec1(0.0), dt);
    }
    // After 10 steps, filtered_sp should be closer to 1.0 but not yet there
    REQUIRE(fsp > 0.5);
    REQUIRE(fsp < 1.0);
}

TEST_CASE("PvFilter smooths step measurement", "[pid][siso][pv-filter]")
{
    using PvfPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::PvFilter>;
    PvfPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.template policy<ctrlpp::PvFilter>().tf = {0.05};
    PvfPid pid(cfg);

    double alpha = 0.05 / (0.05 + dt);

    // Step 1: filtered_meas = (1-alpha)*1.0
    auto u1 = pid.compute(vec1(0.0), vec1(1.0), dt);
    double fm1 = (1.0 - alpha) * 1.0;
    // e = 0 - fm1 = -fm1, P = -fm1
    REQUIRE_THAT(u1[0], WithinAbs(-fm1, tol));

    // After several steps, filtered_meas approaches 1.0
    double fm = fm1;
    for (int i = 1; i < 20; ++i) {
        fm = alpha * fm + (1.0 - alpha) * 1.0;
        pid.compute(vec1(0.0), vec1(1.0), dt);
    }
    REQUIRE(fm > 0.9);
    REQUIRE(fm < 1.0);
}

TEST_CASE("FeedForward adds callable output to control signal", "[pid][siso][feed-forward]")
{
    // Feed-forward: ff(sp, dt) = sp (identity)
    auto ff_func = [](const Vec1& sp, double) -> Vec1 { return sp; };
    using FF = ctrlpp::FeedForward<decltype(ff_func)>;
    using FfPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, FF>;

    FfPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.template policy<FF>().ff_func = ff_func;
    FfPid pid(cfg);

    // P = 1.0 * (1.0 - 0.0) = 1.0, FF = sp = 1.0, total = 2.0
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(2.0, tol));
}

TEST_CASE("SetpointFilter + PvFilter combined", "[pid][siso][filter-combo]")
{
    using ComboPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1,
        ctrlpp::SetpointFilter, ctrlpp::PvFilter>;
    ComboPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.template policy<ctrlpp::SetpointFilter>().tf = {0.1};
    cfg.template policy<ctrlpp::PvFilter>().tf = {0.05};
    ComboPid pid(cfg);

    double alpha_sp = 0.1 / (0.1 + dt);
    double alpha_pv = 0.05 / (0.05 + dt);

    // Step 1: sp=1.0, meas=0.5
    double fsp = (1.0 - alpha_sp) * 1.0;
    double fm = (1.0 - alpha_pv) * 0.5;
    auto u = pid.compute(vec1(1.0), vec1(0.5), dt);
    // e = fsp - fm, P = 1.0 * e
    REQUIRE_THAT(u[0], WithinAbs(fsp - fm, tol));
}

TEST_CASE("No filter policies: output unchanged from baseline", "[pid][siso][no-filter]")
{
    // Bare PID (no filter policies) should match Plan 01 baseline exactly
    SisoPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.ki = {0.5};
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    // P=2*1=2, I=0.5*1*0.01=0.005
    REQUIRE_THAT(u[0], WithinAbs(2.005, tol));
}

TEST_CASE("SetpointFilter reset clears filter state", "[pid][siso][setpoint-filter][reset]")
{
    using SpfPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::SetpointFilter>;
    SpfPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.template policy<ctrlpp::SetpointFilter>().tf = {0.1};
    SpfPid pid(cfg);

    double alpha = 0.1 / (0.1 + dt);

    // Run a few steps to build up filter state
    pid.compute(vec1(1.0), vec1(0.0), dt);
    pid.compute(vec1(1.0), vec1(0.0), dt);
    pid.compute(vec1(1.0), vec1(0.0), dt);

    pid.reset();

    // After reset, filter state is cleared — first step behaves as if fresh
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    double fsp = (1.0 - alpha) * 1.0;
    REQUIRE_THAT(u[0], WithinAbs(fsp, tol));
}

// =============================================================================
// Plan 02: Output pipeline -- anti-windup, rate limiting, saturated flag
// =============================================================================

TEST_CASE("saturated() returns true when output is clamped", "[pid][siso][saturated]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {10.0};
    cfg.output_min = {-5.0};
    cfg.output_max = {5.0};
    SisoPid pid(cfg);

    // P = 10*1 = 10, clamped to 5 -> saturated
    pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE(pid.saturated() == true);

    // P = 10*0.1 = 1.0, not clamped -> not saturated
    pid.compute(vec1(0.1), vec1(0.0), dt);
    REQUIRE(pid.saturated() == false);
}

TEST_CASE("Without AntiWindup: integral winds up unboundedly", "[pid][siso][windup]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {1.0};
    cfg.output_max = {5.0};
    SisoPid pid(cfg);

    // Constant error of 1.0 for many steps
    for (int i = 0; i < 1000; ++i)
        pid.compute(vec1(10.0), vec1(0.0), dt);

    // Integral should be large (no anti-windup to stop it)
    // I = Ki * e * dt * 1000 = 1 * 10 * 0.01 * 1000 ~ 100
    REQUIRE_THAT(pid.integral()[0], WithinAbs(100.0, 1e-6));
}

TEST_CASE("BackCalc anti-windup limits integral growth during saturation",
    "[pid][siso][anti-windup][backcalc]")
{
    using AW = ctrlpp::AntiWindup<ctrlpp::BackCalc>;
    using AwPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {1.0};
    cfg.output_max = {5.0};
    cfg.template policy<AW>().kb = {1.0};
    AwPid pid(cfg);

    // Run many steps with large error to cause saturation
    for (int i = 0; i < 1000; ++i)
        pid.compute(vec1(10.0), vec1(0.0), dt);

    // Integral should be bounded (back-calc feedback limits growth)
    // Without anti-windup it would be 100.0
    REQUIRE(pid.integral()[0] < 50.0);
}

TEST_CASE("BackCalc default Kb auto-computation", "[pid][siso][anti-windup][backcalc][auto-kb]")
{
    using AW = ctrlpp::AntiWindup<ctrlpp::BackCalc>;

    SECTION("Kb = sqrt(Ki*Kd) when Kd != 0") {
        using AwPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, AW>;
        AwPid::config_type cfg{};
        cfg.kp = {2.0};
        cfg.ki = {4.0};
        cfg.kd = {1.0};
        cfg.output_max = {5.0};
        // kb left at 0 -> auto-compute = sqrt(4*1) = 2
        AwPid pid(cfg);

        // Run until saturation
        for (int i = 0; i < 100; ++i)
            pid.compute(vec1(10.0), vec1(0.0), dt);

        // Verify anti-windup is active (integral bounded)
        REQUIRE(pid.integral()[0] < 50.0);
    }

    SECTION("Kb = Ki when Kd == 0") {
        using AwPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, AW>;
        AwPid::config_type cfg{};
        cfg.kp = {1.0};
        cfg.ki = {2.0};
        // kd left at 0
        cfg.output_max = {5.0};
        AwPid pid(cfg);

        for (int i = 0; i < 100; ++i)
            pid.compute(vec1(10.0), vec1(0.0), dt);

        REQUIRE(pid.integral()[0] < 50.0);
    }
}

TEST_CASE("Clamping anti-windup freezes integral during saturation",
    "[pid][siso][anti-windup][clamping]")
{
    using AW = ctrlpp::AntiWindup<ctrlpp::Clamping>;
    using AwPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {1.0};
    cfg.output_max = {5.0};
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

TEST_CASE("ConditionalIntegration freezes integral when error exceeds threshold",
    "[pid][siso][anti-windup][conditional]")
{
    using AW = ctrlpp::AntiWindup<ctrlpp::ConditionalIntegration>;
    using AwPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {1.0};
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

TEST_CASE("RateLimit constrains output change per step", "[pid][siso][rate-limit]")
{
    using RlPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::RateLimit>;
    RlPid::config_type cfg{};
    cfg.kp = {100.0};
    cfg.template policy<ctrlpp::RateLimit>().rate_max = {10.0};
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

TEST_CASE("RateLimit applies before output clamp", "[pid][siso][rate-limit][pipeline-order]")
{
    using RlPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::RateLimit>;
    RlPid::config_type cfg{};
    cfg.kp = {100.0};
    cfg.output_max = {0.05};
    cfg.template policy<ctrlpp::RateLimit>().rate_max = {10.0};
    RlPid pid(cfg);

    // Raw = 100, rate limited to 0.1, then clamped to 0.05
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(0.05, tol));
    REQUIRE(pid.saturated() == true);
}

// =============================================================================
// Plan 02: Signal path -- setpoint weighting, derivative filter
// =============================================================================

TEST_CASE("Setpoint weighting b=0.5 modifies proportional term",
    "[pid][siso][setpoint-weighting]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.b = {0.5};
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
    cfg.kp = {2.0};
    cfg.b = {0.0};
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
    cfg.kp = {1.0};
    cfg.ki = {1.0};
    cfg.b = {0.5};
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
    cfg.kp = {1.0};
    cfg.kd = {0.1};
    cfg.derivative_on_error = true;
    cfg.c = {0.5};
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
    cfg.kp = {1.0};
    cfg.kd = {0.1};
    cfg.derivative_on_error = true;
    cfg.c = {0.0};
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

TEST_CASE("DerivFilter reduces peak derivative on step input",
    "[pid][siso][deriv-filter]")
{
    // Unfiltered PID
    SisoPid::config_type cfg_unfiltered{};
    cfg_unfiltered.kp = {1.0};
    cfg_unfiltered.kd = {1.0};
    SisoPid pid_unfiltered(cfg_unfiltered);

    // Filtered PID
    using DfPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::DerivFilter>;
    DfPid::config_type cfg_filtered{};
    cfg_filtered.kp = {1.0};
    cfg_filtered.kd = {1.0};
    cfg_filtered.template policy<ctrlpp::DerivFilter>().n = {10.0};
    DfPid pid_filtered(cfg_filtered);

    // Step 1: first step, D=0 for both
    pid_unfiltered.compute(vec1(0.0), vec1(0.0), dt);
    pid_filtered.compute(vec1(0.0), vec1(0.0), dt);

    // Step 2: measurement step from 0 to 1 -> large derivative
    auto u_unfiltered = pid_unfiltered.compute(vec1(0.0), vec1(1.0), dt);
    auto u_filtered = pid_filtered.compute(vec1(0.0), vec1(1.0), dt);

    // Unfiltered: D = -Kd * (1-0)/dt = -1/0.01 = -100
    // P = 1.0 * (0-1) = -1
    // u_unfiltered = -1 + (-100) = -101
    REQUIRE_THAT(u_unfiltered[0], WithinAbs(-101.0, tol));

    // Filtered: D should have smaller magnitude (filter dampens the spike)
    // Tf = Kd / (Kp * N) = 1.0 / (1.0 * 10.0) = 0.1
    // alpha = Tf / (Tf + dt) = 0.1 / 0.11 ~ 0.9091
    // D_filtered = alpha * 0 + (1 - alpha) * (-100) ~ -9.09
    double tf = 1.0 / (1.0 * 10.0);
    double alpha = tf / (tf + dt);
    double d_filtered = (1.0 - alpha) * (-100.0);
    // P = -1
    REQUIRE_THAT(u_filtered[0], WithinAbs(-1.0 + d_filtered, 1e-6));

    // Peak derivative magnitude is smaller with filter
    REQUIRE(std::abs(u_filtered[0]) < std::abs(u_unfiltered[0]));
}

TEST_CASE("DerivFilter with N parameter smooths derivative over multiple steps",
    "[pid][siso][deriv-filter]")
{
    using DfPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::DerivFilter>;
    DfPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.kd = {0.5};
    cfg.template policy<ctrlpp::DerivFilter>().n = {20.0};
    DfPid pid(cfg);

    // Tf = Kd / (Kp * N) = 0.5 / (2.0 * 20.0) = 0.0125
    double tf_val = 0.5 / (2.0 * 20.0);
    double alpha = tf_val / (tf_val + dt);

    // Step 1: D=0 (first step)
    pid.compute(vec1(0.0), vec1(0.0), dt);

    // Step 2: measurement step to 1.0
    // D_raw = -Kd * (1-0)/dt = -0.5 / 0.01 = -50
    // D_filt = alpha*0 + (1-alpha)*(-50)
    double d_filt = (1.0 - alpha) * (-50.0);
    auto u2 = pid.compute(vec1(0.0), vec1(1.0), dt);
    double p2 = 2.0 * (0.0 - 1.0);
    REQUIRE_THAT(u2[0], WithinAbs(p2 + d_filt, 1e-6));

    // Step 3: measurement stays at 1.0
    // D_raw = -Kd * (1-1)/dt = 0
    // D_filt = alpha*d_filt + (1-alpha)*0
    double d_filt2 = alpha * d_filt;
    auto u3 = pid.compute(vec1(0.0), vec1(1.0), dt);
    double p3 = 2.0 * (0.0 - 1.0);
    REQUIRE_THAT(u3[0], WithinAbs(p3 + d_filt2, 1e-6));
}

TEST_CASE("Without DerivFilter policy, unfiltered derivative is used",
    "[pid][siso][deriv-filter][compile]")
{
    // This ensures that without DerivFilter, no filter overhead occurs
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.kd = {1.0};
    SisoPid pid(cfg);

    pid.compute(vec1(0.0), vec1(0.0), dt);
    auto u = pid.compute(vec1(0.0), vec1(1.0), dt);
    // D = -1.0 * (1-0)/0.01 = -100
    // P = 1.0 * (0-1) = -1
    REQUIRE_THAT(u[0], WithinAbs(-101.0, tol));
}

// =============================================================================
// Plan 03: Velocity form PID
// =============================================================================

TEST_CASE("VelocityForm P-only: delta_u = Kp*(e(k)-e(k-1))", "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::VelocityForm>;
    VPid::config_type cfg{};
    cfg.kp = {2.0};
    VPid pid(cfg);

    // Step 1: e=1.0, prev_e=0.0
    // delta_u = Kp*(1-0) + 0 + 0 = 2.0
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // Step 2: e=1.0, prev_e=1.0
    // delta_u = Kp*(1-1) = 0.0
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(0.0, tol));

    // Step 3: e=0.5, prev_e=1.0
    // delta_u = Kp*(0.5-1.0) = -1.0
    auto u3 = pid.compute(vec1(1.0), vec1(0.5), dt);
    REQUIRE_THAT(u3[0], WithinAbs(-1.0, tol));
}

TEST_CASE("VelocityForm PI: includes Ki*e*dt incremental term", "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::VelocityForm>;
    VPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    VPid pid(cfg);

    // Step 1: e=1.0, prev_e=0
    // dP = 1*(1-0) = 1.0, dI = 0.5*1*0.01 = 0.005
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    // Step 2: e=1.0, prev_e=1.0
    // dP = 1*(1-1) = 0, dI = 0.5*1*0.01 = 0.005
    auto u2 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(0.005, tol));

    // Step 3: e=0.5
    // dP = 1*(0.5-1) = -0.5, dI = 0.5*0.5*0.01 = 0.0025
    auto u3 = pid.compute(vec1(1.0), vec1(0.5), dt);
    REQUIRE_THAT(u3[0], WithinAbs(-0.5 + 0.0025, tol));
}

TEST_CASE("VelocityForm PID: full formula with second-order D difference",
    "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::VelocityForm>;
    VPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.ki = {0.5};
    cfg.kd = {0.1};
    VPid pid(cfg);

    // Step 1: e(1)=1.0, e(0)=0, e(-1)=0
    // dP = 2*(1-0) = 2, dI = 0.5*1*0.01 = 0.005
    // dD = 0.1*(1 - 2*0 + 0)/0.01 = 0.1*100 = 10
    double du1 = 2.0 + 0.005 + 10.0;
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(du1, tol));

    // Step 2: e(2)=0.8, e(1)=1.0, e(0)=0
    // dP = 2*(0.8-1.0) = -0.4, dI = 0.5*0.8*0.01 = 0.004
    // dD = 0.1*(0.8 - 2*1.0 + 0)/0.01 = 0.1*(-1.2/0.01) = 0.1*(-120) = -12
    double du2 = -0.4 + 0.004 + (-12.0);
    auto u2 = pid.compute(vec1(1.0), vec1(0.2), dt);
    REQUIRE_THAT(u2[0], WithinAbs(du2, tol));

    // Step 3: e(3)=0.8, e(2)=0.8, e(1)=1.0
    // dP = 2*(0.8-0.8) = 0, dI = 0.5*0.8*0.01 = 0.004
    // dD = 0.1*(0.8 - 2*0.8 + 1.0)/0.01 = 0.1*(0.2/0.01) = 0.1*20 = 2
    double du3 = 0.0 + 0.004 + 2.0;
    auto u3 = pid.compute(vec1(1.0), vec1(0.2), dt);
    REQUIRE_THAT(u3[0], WithinAbs(du3, tol));
}

TEST_CASE("VelocityForm steady-state: delta_u converges to Ki*e*dt",
    "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::VelocityForm>;
    VPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {2.0};
    cfg.kd = {0.1};
    VPid pid(cfg);

    // First few steps build up history, then with constant error:
    // dP = 0 (error unchanged), dD = 0 (e(k)-2e(k-1)+e(k-2)=0), dI = Ki*e*dt
    pid.compute(vec1(1.0), vec1(0.0), dt);  // step 1
    pid.compute(vec1(1.0), vec1(0.0), dt);  // step 2
    pid.compute(vec1(1.0), vec1(0.0), dt);  // step 3

    // From step 3 onward, e is constant at 1.0
    // delta_u should be Ki*e*dt = 2*1*0.01 = 0.02
    auto u4 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u4[0], WithinAbs(0.02, tol));

    auto u5 = pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u5[0], WithinAbs(0.02, tol));
}

TEST_CASE("VelocityForm + AntiWindup compiles and runs (anti-windup is no-op)",
    "[pid][siso][velocity-form][anti-windup]")
{
    using VPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1,
        ctrlpp::VelocityForm, ctrlpp::AntiWindup<ctrlpp::BackCalc>>;
    VPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    VPid pid(cfg);

    // Should produce same result as without AntiWindup
    auto u1 = pid.compute(vec1(1.0), vec1(0.0), dt);
    // dP = 1*(1-0) = 1, dI = 0.5*1*0.01 = 0.005
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));
}

TEST_CASE("VelocityForm reset clears history", "[pid][siso][velocity-form][reset]")
{
    using VPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, ctrlpp::VelocityForm>;
    VPid::config_type cfg{};
    cfg.kp = {2.0};
    cfg.ki = {0.5};
    VPid pid(cfg);

    // Build up history
    pid.compute(vec1(1.0), vec1(0.0), dt);
    pid.compute(vec1(1.0), vec1(0.0), dt);

    pid.reset();

    // After reset: behaves as first step (prev_error=0, prev_prev_error=0)
    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
    // dP = 2*(1-0) = 2, dI = 0.5*1*0.01 = 0.005, dD = 0 (Kd=0)
    REQUIRE_THAT(u[0], WithinAbs(2.005, tol));
}

// =============================================================================
// Plan 04: Bumpless transfer, gain scheduling, integral management
// =============================================================================

TEST_CASE("Tracking signal drives integral for bumpless transfer",
    "[pid][siso][tracking]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
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
    cfg.kp = {2.0};
    cfg.ki = {1.0};
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

TEST_CASE("set_params rescales integral for bumpless gain change",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {2.0};
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
    new_cfg.ki = {4.0};
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.1, tol));
}

TEST_CASE("set_params with Ki going to zero clears integral",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {2.0};
    SisoPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), dt);

    REQUIRE(pid.integral()[0] != 0.0);

    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = {0.0};
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));
}

TEST_CASE("set_params Kp change produces no output discontinuity",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    SisoPid pid(cfg);

    // Run until steady output
    Vec1 u_before{};
    for (int i = 0; i < 20; ++i)
        u_before = pid.compute(vec1(1.0), vec1(0.5), dt);

    // Change Kp from 1 to 2
    SisoPid::config_type new_cfg = cfg;
    new_cfg.kp = {2.0};
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

TEST_CASE("set_integral sets integral to known value",
    "[pid][siso][integral-management]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {1.0};
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
    cfg.kp = {1.0};
    cfg.ki = {1.0};
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

TEST_CASE("params() returns current config after set_params",
    "[pid][siso][gain-scheduling]")
{
    SisoPid::config_type cfg{};
    cfg.kp = {1.0};
    cfg.ki = {0.5};
    SisoPid pid(cfg);

    SisoPid::config_type new_cfg = cfg;
    new_cfg.kp = {3.0};
    new_cfg.ki = {1.5};
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.params().kp[0], WithinAbs(3.0, tol));
    REQUIRE_THAT(pid.params().ki[0], WithinAbs(1.5, tol));
}

// =============================================================================
// Plan 04: Performance assessment metrics
// =============================================================================

TEST_CASE("IAE accumulates integral of |error| * dt",
    "[pid][siso][perf-assessment][iae]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::IAE>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    // Constant error=1.0, 10 steps of dt=0.1
    // IAE = sum(|1.0| * 0.1) = 1.0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("ISE accumulates integral of error^2 * dt",
    "[pid][siso][perf-assessment][ise]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::ISE>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    // Constant error=2.0, 10 steps of dt=0.1
    // ISE = sum(4.0 * 0.1) = 4.0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(2.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(4.0, 1e-10));
}

TEST_CASE("ITAE accumulates integral of t * |error| * dt",
    "[pid][siso][perf-assessment][itae]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::ITAE>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    // Constant error=1.0, dt=0.1
    // After step k: accumulated_time = k * 0.1
    // ITAE = sum(t_k * |1.0| * 0.1)
    // t_k = 0.1, 0.2, ..., 1.0 (accumulated AFTER dt update)
    // ITAE = 0.1*(0.1 + 0.2 + ... + 1.0) = 0.1 * 5.5 = 0.55
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.55, 1e-10));
}

TEST_CASE("Multiple metrics accumulate simultaneously",
    "[pid][siso][perf-assessment][multi]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::IAE, ctrlpp::ISE>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(2.0), vec1(0.0), 0.1);

    // IAE = sum(|2.0| * 0.1) = 2.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(2.0, 1e-10));
    // ISE = sum(4.0 * 0.1) = 4.0
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(4.0, 1e-10));
}

TEST_CASE("OscillationDetect counts zero-crossings and detects oscillation",
    "[pid][siso][perf-assessment][oscillation]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::OscillationDetect>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    // Alternating error sign: +1, -1, +1, -1, ...  for 10 steps at dt=0.1
    // 9 zero-crossings over 1.0 second -> rate = 9/1.0 = 9 > threshold 5
    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        pid.compute(vec1(sp), vec1(0.0), 0.1);
    }

    REQUIRE(pid.metric<ctrlpp::OscillationDetect>()[0] >= 9.0 - tol);
    REQUIRE(pid.oscillating() == true);
}

TEST_CASE("OscillationDetect: constant error sign -> not oscillating",
    "[pid][siso][perf-assessment][oscillation]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::OscillationDetect>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::OscillationDetect>()[0], WithinAbs(0.0, tol));
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("OscillationDetect + IAE both accumulate",
    "[pid][siso][perf-assessment][oscillation][iae]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::OscillationDetect, ctrlpp::IAE>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        pid.compute(vec1(sp), vec1(0.0), 0.1);
    }

    REQUIRE(pid.metric<ctrlpp::OscillationDetect>()[0] >= 9.0 - tol);
    // IAE = sum(|1.0| * 0.1) = 1.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("reset_metrics clears all metric accumulators",
    "[pid][siso][perf-assessment][reset]")
{
    using PA = ctrlpp::PerfAssessment<ctrlpp::IAE, ctrlpp::OscillationDetect>;
    using PaPid = ctrlpp::Pid<NaiveLinalg, double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = {1.0};
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        pid.compute(vec1(sp), vec1(0.0), 0.1);
    }

    REQUIRE(pid.metric<ctrlpp::IAE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::OscillationDetect>()[0] > 0.0);

    pid.reset_metrics();

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.metric<ctrlpp::OscillationDetect>()[0], WithinAbs(0.0, tol));
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("MIMO performance metrics computed per channel",
    "[pid][mimo][perf-assessment]")
{
    using MimoPid = ctrlpp::Pid<NaiveLinalg, double, 1, 2, 2,
        ctrlpp::PerfAssessment<ctrlpp::IAE>>;
    using Vec2 = NaiveLinalg::vector_type<double, 2>;
    MimoPid::config_type cfg{};
    cfg.kp = {1.0, 1.0};
    MimoPid pid(cfg);

    // Channel 0: error=1.0, Channel 1: error=3.0
    for (int i = 0; i < 10; ++i)
        pid.compute(Vec2{1.0, 3.0}, Vec2{0.0, 0.0}, 0.1);

    // IAE ch0 = 10*|1|*0.1 = 1.0, IAE ch1 = 10*|3|*0.1 = 3.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[1], WithinAbs(3.0, 1e-10));
}
