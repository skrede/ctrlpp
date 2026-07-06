#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1, 1, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

// Bound for comparing two outputs computed with a different operation order. Each of
// the arithmetic operations composing the output is rounded once per controller, so
// the accumulated forward error is bounded by twice the operation count times one
// unit in the last place at the working scale.
double output_continuity_tol(double scale)
{
    constexpr int arithmetic_ops = 8; // products/sums forming the PID output
    return 2.0 * arithmetic_ops * std::numeric_limits<double>::epsilon() * std::abs(scale);
}

}

// The integral state is stored in output units (each increment is ki*e*dt and enters
// the output directly), so its contribution to the output is already continuous across
// a gain change. set_params must therefore leave the integral state untouched: this is
// what makes gain scheduling bumpless. Rescaling the state by ki_old/ki_new would
// instead step the integral contribution, injecting a bump.
TEST_CASE("set_params with a Ki change is bumpless (integral state unchanged)",
    "[pid][siso][gain-scheduling]")
{
    const double kp = 1.0, ki_old = 2.0, ki_new = 4.0;
    const double sp = 1.0, meas = 0.0;
    const double e = sp - meas; // setpoint weight b defaults to 1

    SisoPid::config_type cfg{};
    cfg.kp = vec1(kp);
    cfg.ki = vec1(ki_old);

    // Two controllers driven identically hold bit-identical state; one then changes Ki.
    SisoPid ref(cfg);
    SisoPid sched(cfg);
    for (int i = 0; i < 10; ++i) {
        ref.compute(vec1(sp), vec1(meas), Ts);
        sched.compute(vec1(sp), vec1(meas), Ts);
    }

    const double integral_before = sched.integral()[0];
    REQUIRE_THAT(integral_before, WithinAbs(ref.integral()[0], tol));

    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(ki_new);
    sched.set_params(new_cfg);

    // Bumpless transfer: the stored integral (in output units) is not rescaled.
    REQUIRE_THAT(sched.integral()[0], WithinAbs(integral_before, tol));

    // On the next identical step the outputs differ only by the intended change in the
    // integral increment on the current error, (ki_new - ki_old)*e*dt, with no jump in
    // the accumulated integral contribution.
    const double u_ref = ref.compute(vec1(sp), vec1(meas), Ts)[0];
    const double u_sched = sched.compute(vec1(sp), vec1(meas), Ts)[0];
    const double expected_diff = (ki_new - ki_old) * e * Ts;

    REQUIRE_THAT(u_sched - u_ref, WithinAbs(expected_diff, output_continuity_tol(u_sched)));
}

TEST_CASE("set_params with Ki going to zero preserves the accumulated integral",
    "[pid][siso][gain-scheduling]")
{
    const double kp = 1.0, ki = 2.0;
    const double sp = 1.0, meas = 0.0;

    SisoPid::config_type cfg{};
    cfg.kp = vec1(kp);
    cfg.ki = vec1(ki);
    SisoPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(sp), vec1(meas), Ts);

    const double integral_before = pid.integral()[0];
    REQUIRE(integral_before != 0.0);

    // Turning integral action off must not clear the accumulated integral: that would
    // step the output. The integrator simply stops growing (increment is now zero)
    // while its stored contribution persists.
    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(0.0);
    pid.set_params(new_cfg);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_before, tol));

    // A further step adds no increment and holds the integral, so the output is the
    // proportional term plus the retained integral contribution.
    const double u = pid.compute(vec1(sp), vec1(meas), Ts)[0];
    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_before, tol));
    REQUIRE_THAT(u, WithinAbs(kp * (sp - meas) + integral_before, output_continuity_tol(u)));
}

TEST_CASE("set_params with a Kp change steps only the proportional term",
    "[pid][siso][gain-scheduling]")
{
    const double kp_old = 1.0, kp_new = 2.0, ki = 0.5;
    const double sp = 1.0, meas = 0.5;
    const double ep = sp - meas; // setpoint weight b defaults to 1

    SisoPid::config_type cfg{};
    cfg.kp = vec1(kp_old);
    cfg.ki = vec1(ki);

    SisoPid ref(cfg);
    SisoPid sched(cfg);
    for (int i = 0; i < 20; ++i) {
        ref.compute(vec1(sp), vec1(meas), Ts);
        sched.compute(vec1(sp), vec1(meas), Ts);
    }

    const double integral_before = sched.integral()[0];

    SisoPid::config_type new_cfg = cfg;
    new_cfg.kp = vec1(kp_new);
    sched.set_params(new_cfg);

    // Changing Kp leaves the integral state untouched (Ki unchanged), so the only
    // difference on the next step is the deliberate proportional change on the current
    // error, (kp_new - kp_old)*ep.
    REQUIRE_THAT(sched.integral()[0], WithinAbs(integral_before, tol));

    const double u_ref = ref.compute(vec1(sp), vec1(meas), Ts)[0];
    const double u_sched = sched.compute(vec1(sp), vec1(meas), Ts)[0];
    const double expected_diff = (kp_new - kp_old) * ep;

    REQUIRE_THAT(u_sched - u_ref, WithinAbs(expected_diff, output_continuity_tol(u_sched)));
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
    auto u = pid.compute(vec1(1.0), vec1(0.0), Ts);
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
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    double integral_val = pid.integral()[0];

    // Freeze
    pid.freeze_integral(true);

    // Run 10 more steps -- integral should not change
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), Ts);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_val, tol));

    // Unfreeze
    pid.freeze_integral(false);

    // Integral should resume
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE(pid.integral()[0] > integral_val);
}

TEST_CASE("set_params from Ki=0 to Ki!=0 leaves the zero integral in place",
    "[pid][siso][gain-scheduling][edge-case]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.0);
    SisoPid pid(cfg);

    // Run some steps with Ki=0 -> integral never accumulates, stays 0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Enabling integral action does not touch the (zero) integral state; it simply
    // starts accumulating from where it was.
    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(2.0);
    pid.set_params(new_cfg);

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Now integral should accumulate with new Ki
    pid.compute(vec1(1.0), vec1(0.0), Ts);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(2.0 * 1.0 * Ts, tol));
}

TEST_CASE("ISA form set_params with a Ti change is bumpless",
    "[pid][siso][gain-scheduling][isa]")
{
    using IsaPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::isa_form>;
    const double kp = 2.0, ti_old = 4.0, ti_new = 2.0;
    const double sp = 1.0, meas = 0.0;
    const double e = sp - meas;
    // Internal integral gain in ISA form is Kp/Ti.
    const double ki_int_old = kp / ti_old; // 0.5
    const double ki_int_new = kp / ti_new; // 1.0

    IsaPid::config_type cfg{};
    cfg.kp = vec1(kp);
    cfg.ki = vec1(ti_old);
    cfg.kd = vec1(0.0);

    IsaPid ref(cfg);
    IsaPid sched(cfg);
    for (int i = 0; i < 10; ++i) {
        ref.compute(vec1(sp), vec1(meas), Ts);
        sched.compute(vec1(sp), vec1(meas), Ts);
    }

    const double integral_before = sched.integral()[0];
    REQUIRE_THAT(integral_before, WithinAbs(ref.integral()[0], tol));

    IsaPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(ti_new);
    sched.set_params(new_cfg);

    // The stored integral (output units) is not rescaled when Ti changes.
    REQUIRE_THAT(sched.integral()[0], WithinAbs(integral_before, tol));

    const double u_ref = ref.compute(vec1(sp), vec1(meas), Ts)[0];
    const double u_sched = sched.compute(vec1(sp), vec1(meas), Ts)[0];
    const double expected_diff = (ki_int_new - ki_int_old) * e * Ts;

    REQUIRE_THAT(u_sched - u_ref, WithinAbs(expected_diff, output_continuity_tol(u_sched)));
}

TEST_CASE("repeated set_params never rescales the integral state",
    "[pid][siso][gain-scheduling][edge-case]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    // Accumulate integral = Ki*e*dt*5 = 1.0*1.0*0.01*5 = 0.05
    for (int i = 0; i < 5; ++i)
        pid.compute(vec1(1.0), vec1(0.0), Ts);
    const double integral_before = pid.integral()[0];
    REQUIRE_THAT(integral_before, WithinAbs(0.05, tol));

    // Doubling Ki does not halve the stored integral: the state is in output units and
    // is left untouched, so the integral contribution stays continuous.
    SisoPid::config_type new_cfg = cfg;
    new_cfg.ki = vec1(2.0);
    pid.set_params(new_cfg);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_before, tol));

    // Halving Ki back likewise leaves the stored integral in place.
    SisoPid::config_type newer_cfg = cfg;
    newer_cfg.ki = vec1(1.0);
    pid.set_params(newer_cfg);
    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_before, tol));
}
