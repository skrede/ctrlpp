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

TEST_CASE("P-only controller", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.5);
    SisoPid pid(cfg);

    auto u = pid.compute(vec1(1.0), vec1(0.0), dt);
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
    using IsaPid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::isa_form>;
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

    auto u_isa = isa_pid.compute(vec1(1.0), vec1(0.0), dt);
    auto u_par = par_pid.compute(vec1(1.0), vec1(0.0), dt);
    REQUIRE_THAT(u_isa[0], WithinAbs(u_par[0], tol));

    auto u_isa2 = isa_pid.compute(vec1(1.0), vec1(0.1), dt);
    auto u_par2 = par_pid.compute(vec1(1.0), vec1(0.1), dt);
    REQUIRE_THAT(u_isa2[0], WithinAbs(u_par2[0], tol));
}

TEST_CASE("reset() zeros state", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(1.0);
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
    cfg.kp = vec1(0.0);
    cfg.kd = vec1(10.0);
    SisoPid pid(cfg);

    // Even with large measurement, derivative should be zero on first step
    auto u = pid.compute(vec1(0.0), vec1(100.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(0.0, tol));
}

TEST_CASE("Zero dt returns previous output", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
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

TEST_CASE("error() returns last error", "[pid][siso]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    pid.compute(vec1(3.0), vec1(1.0), dt);
    REQUIRE_THAT(pid.error()[0], WithinAbs(2.0, tol));
}
