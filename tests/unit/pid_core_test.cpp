#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

}

TEST_CASE("P-only controller", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.5);
    SisoPid pid(cfg);

    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u[0], WithinAbs(2.5, tol));
}

TEST_CASE("Full PID matches hand computation", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    cfg.kd = vec1(0.1);
    SisoPid pid(cfg);

    // Step 1: sp=1, meas=0, e=1
    // P = 2*1 = 2, I = 0.5*1*0.01 = 0.005, D = 0 (first step)
    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(2.005, tol));

    // Step 2: sp=1, meas=0.1, e=0.9
    // P = 2*0.9 = 1.8
    // I = 0.005 + 0.5*0.9*0.01 = 0.0095
    // D = -0.1*(0.1-0.0)/0.01 = -1.0 (derivative on measurement)
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.1), Ts));
    REQUIRE_THAT(u2[0], WithinAbs(1.8 + 0.0095 - 1.0, tol));
}

TEST_CASE("ISA form converts to parallel form", "[pid][siso][isa]")
{
    // ISA form: Kp=2, Ti=4, Td=0.05
    // Parallel: Kp=2, Ki=2/4=0.5, Kd=2*0.05=0.1
    using IsaPid = ctrlpp::pid<double, 1, ctrlpp::isa_form>;
    IsaPid::config_type isa_cfg{};
    isa_cfg.kp = vec1(2.0);
    isa_cfg.ki = vec1(4.0);  // This is Ti for ISA form
    isa_cfg.kd = vec1(0.05); // This is Td for ISA form
    IsaPid isa_pid(isa_cfg);

    // Parallel form reference
    SisoPid::config_type par_cfg{};
    par_cfg.kp = vec1(2.0);
    par_cfg.ki = vec1(0.5);
    par_cfg.kd = vec1(0.1);
    SisoPid par_pid(par_cfg);

    auto u_isa = ctrlpp::test::commanded(isa_pid.compute(vec1(1.0), vec1(0.0), Ts));
    auto u_par = ctrlpp::test::commanded(par_pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u_isa[0], WithinAbs(u_par[0], tol));

    auto u_isa2 = ctrlpp::test::commanded(isa_pid.compute(vec1(1.0), vec1(0.1), Ts));
    auto u_par2 = ctrlpp::test::commanded(par_pid.compute(vec1(1.0), vec1(0.1), Ts));
    REQUIRE_THAT(u_isa2[0], WithinAbs(u_par2[0], tol));
}

TEST_CASE("reset() zeros state", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());
    REQUIRE(pid.integral()[0] != 0.0);

    pid.reset();
    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.error()[0], WithinAbs(0.0, tol));

    // After reset, first step derivative should be zero again
    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.5), Ts));
    // P = 2*0.5 = 1.0, I = 1*0.5*0.01 = 0.005, D = 0 (first step after reset)
    REQUIRE_THAT(u[0], WithinAbs(1.005, tol));
}

TEST_CASE("First step derivative is zero", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.0);
    cfg.kd = vec1(10.0);
    SisoPid pid(cfg);

    // Even with large measurement, derivative should be zero on first step
    auto u = ctrlpp::test::commanded(pid.compute(vec1(0.0), vec1(100.0), Ts));
    REQUIRE_THAT(u[0], WithinAbs(0.0, tol));
}

TEST_CASE("Zero and negative dt are rejected rather than answered with the held output",
    "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    SisoPid pid(cfg);

    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // A stopped clock used to be answered with the stored output on the success
    // path, which a caller could not distinguish from a command the controller
    // had actually computed. It now names the fault instead, and the caller
    // decides what the actuator does with a cycle that produced nothing.
    auto stopped = pid.compute(vec1(5.0), vec1(0.0), 0.0);
    REQUIRE_FALSE(stopped.has_value());
    CHECK(stopped.error() == ctrlpp::pid_step_error::invalid_timestep);

    auto backwards = pid.compute(vec1(5.0), vec1(0.0), -0.01);
    REQUIRE_FALSE(backwards.has_value());
    CHECK(backwards.error() == ctrlpp::pid_step_error::invalid_timestep);

    // Neither rejected cycle touched the carried state, so the controller still
    // answers the next valid cycle exactly as it would have.
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    CHECK(u2 == u1);
    CHECK(pid.health() == ctrlpp::pid_health::ok);
}

TEST_CASE("error() returns last error", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    REQUIRE(pid.compute(vec1(3.0), vec1(1.0), Ts).has_value());
    REQUIRE_THAT(pid.error()[0], WithinAbs(2.0, tol));
}

TEST_CASE("ISA form with Ti=0 sets Ki to zero (no integral action)",
    "[pid][siso][isa][edge-case]")
{
    using IsaPid = ctrlpp::pid<double, 1, ctrlpp::isa_form>;
    IsaPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.0);  // Ti=0 in ISA form -> Ki should be 0 (no divide-by-zero)
    cfg.kd = vec1(0.05);
    IsaPid pid(cfg);

    // With Ki=0, integral should stay at zero after multiple steps
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));
}

TEST_CASE("Tracking signal with zero dt is rejected without modifying the integral",
    "[pid][siso][tracking][edge-case]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    // Run one normal step
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());
    const auto integral_before = pid.integral();

    // The tracking overload delegates, so the delegated cycle's rejection is
    // forwarded with its own cause and the back-assignment into the integrator
    // never runs. Exact comparison: nothing was written, so nothing rounded.
    auto stopped = pid.compute(vec1(1.0), vec1(0.0), 0.0, vec1(10.0));
    REQUIRE_FALSE(stopped.has_value());
    CHECK(stopped.error() == ctrlpp::pid_step_error::invalid_timestep);
    CHECK(pid.integral() == integral_before);

    auto backwards = pid.compute(vec1(1.0), vec1(0.0), -0.01, vec1(10.0));
    REQUIRE_FALSE(backwards.has_value());
    CHECK(backwards.error() == ctrlpp::pid_step_error::invalid_timestep);
    CHECK(pid.integral() == integral_before);

    // A non-finite tracking signal is caught before the cycle is delegated, so
    // it cannot destroy the integrator after the rest of the cycle committed.
    auto poisoned = pid.compute(vec1(1.0), vec1(0.0), Ts,
        vec1(std::numeric_limits<double>::quiet_NaN()));
    REQUIRE_FALSE(poisoned.has_value());
    CHECK(poisoned.error() == ctrlpp::pid_step_error::non_finite_tracking_signal);
    CHECK(pid.integral() == integral_before);
    CHECK(pid.error() == vec1(1.0));
}

TEST_CASE("Tracking signal with negative error adjusts integral correctly",
    "[pid][siso][tracking][negative-error]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    SisoPid pid(cfg);

    // sp=0, meas=1 -> negative error
    REQUIRE(pid.compute(vec1(0.0), vec1(1.0), Ts, vec1(2.0)).has_value());

    // p = kp * (b*sp - meas) = 1*(0 - 1) = -1
    // non_integral = u - integral
    // integral = tracking - non_integral = 2.0 - (u - integral)
    // The integral should be set so that on next auto step, output is near tracking
    // integral after tracking = tracking - p = 2.0 - (-1.0) = 3.0
    REQUIRE_THAT(pid.integral()[0], WithinAbs(3.0, tol));
}
