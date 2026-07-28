#include "hardening_helpers.h"
#include "ctrlpp/control/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

using SisoPid = ctrlpp::pid<double, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

Vec1 vec1(double v)
{
    Vec1 r;
    r << v;
    return r;
}

constexpr double Ts = 0.01;

}

namespace {

// A controller with integral and derivative action, so a rejected cycle has
// carried state worth asserting about. A proportional-only controller carries
// nothing the poison could latch into, which is why asserting on its output
// said nothing.
auto make_pid_with_memory() -> SisoPid
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(2.0);
    cfg.kd = vec1(0.5);
    return SisoPid{cfg};
}

}

TEST_CASE("PID NaN setpoint is rejected without touching the carried state",
          "[pid][hardening][negative]")
{
    auto pid = make_pid_with_memory();
    // Drive one valid cycle first, so the integrator and the error history hold
    // something other than zero and "unchanged" is a real claim.
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    const auto integral_before = pid.integral();
    const auto error_before = pid.error();

    auto rejected = pid.compute(vec1(std::numeric_limits<double>::quiet_NaN()), vec1(0.0), Ts);

    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error() == ctrlpp::pid_step_error::non_finite_setpoint);
    // Exact, never a tolerance: a rejected cycle performs no arithmetic on the
    // carried state at all, so bitwise equality is the contract. A tolerance
    // would admit a cycle that partially ran. Both pieces of carried state are
    // asserted -- checking only the integrator would leave the derivative
    // history free to be poisoned.
    CHECK(pid.integral() == integral_before);
    CHECK(pid.error() == error_before);
    // The rejection describes the argument, not the controller, so it must not
    // latch. A latch here would be a library defect, not a strict assertion.
    CHECK(pid.health() == ctrlpp::pid_health::ok);

    // The claim the rejection alone does not establish: the rejected cycle left
    // no trace, so a following valid cycle produces exactly what it would have
    // produced had the rejected one never been attempted.
    auto reference = make_pid_with_memory();
    REQUIRE(reference.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    auto after = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.25), Ts));
    auto expected = ctrlpp::test::commanded(reference.compute(vec1(1.0), vec1(0.25), Ts));
    CHECK(after == expected);
    CHECK(pid.integral() == reference.integral());
    CHECK(pid.error() == reference.error());
}

TEST_CASE("PID NaN measurement is rejected without touching the carried state",
          "[pid][hardening][negative]")
{
    auto pid = make_pid_with_memory();
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    const auto integral_before = pid.integral();
    const auto error_before = pid.error();

    auto rejected = pid.compute(vec1(1.0), vec1(std::numeric_limits<double>::quiet_NaN()), Ts);

    REQUIRE_FALSE(rejected.has_value());
    // A bad measurement names the sensor, a bad setpoint names the reference
    // generator. The two are separated because the caller repairs them in
    // different places.
    CHECK(rejected.error() == ctrlpp::pid_step_error::non_finite_measurement);
    CHECK(pid.integral() == integral_before);
    CHECK(pid.error() == error_before);
    CHECK(pid.health() == ctrlpp::pid_health::ok);

    auto reference = make_pid_with_memory();
    REQUIRE(reference.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    auto after = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.25), Ts));
    auto expected = ctrlpp::test::commanded(reference.compute(vec1(1.0), vec1(0.25), Ts));
    CHECK(after == expected);
    CHECK(pid.integral() == reference.integral());
    CHECK(pid.error() == reference.error());
}

TEST_CASE("PID Inf gains produce finite or Inf output without crash", "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(std::numeric_limits<double>::infinity());
    SisoPid pid(cfg);

    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    CHECK(!std::isnan(u[0]));
}

TEST_CASE("PID near-zero dt returns previous output", "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    SisoPid pid(cfg);

    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), 1e-15));
    CHECK_THAT(u2[0], WithinAbs(u1[0], 1e-12));
}

TEST_CASE("PID P-only exact output", "[pid][hardening][precision]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u[0], WithinAbs(1.0, 1e-15));
}

TEST_CASE("PID precision with known PID gains", "[pid][hardening][precision]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    cfg.kd = vec1(0.0);
    SisoPid pid(cfg);

    // Step 1: e=1, P=2, I=0.5*1*0.01=0.005, D=0
    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u[0], WithinAbs(2.005, 1e-14));
}

TEST_CASE("PID stays bounded over 10000 steps with stable gains", "[pid][hardening][stability]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.5);
    cfg.ki = vec1(0.1);
    cfg.kd = vec1(0.01);
    SisoPid pid(cfg);

    double plant_x = 0.0;
    bool all_finite = true;

    for(int k = 0; k < 10000; ++k)
    {
        auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(plant_x), Ts));
        plant_x = 0.9 * plant_x + 0.1 * u[0];
        if(!std::isfinite(u[0]) || !std::isfinite(plant_x))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
    CHECK(std::abs(plant_x) < 1e6);
}

TEST_CASE("PID with integral action converges to zero error", "[pid][hardening][convergence]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(1.0);
    SisoPid pid(cfg);

    double plant_x = 0.0;

    for(int k = 0; k < 5000; ++k)
    {
        auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(plant_x), Ts));
        plant_x = 0.9 * plant_x + 0.1 * u[0];
    }

    double error = std::abs(1.0 - plant_x);
    REQUIRE(error < 0.01);
}

TEST_CASE("PID extreme gains produce finite output", "[pid][hardening][robustness]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1e12);
    cfg.ki = vec1(1e12);
    SisoPid pid(cfg);

    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    CHECK(std::isfinite(u[0]));

    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.5), Ts));
    CHECK(std::isfinite(u2[0]));
}

TEST_CASE("PID zero gains produce zero output", "[pid][hardening][robustness]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.0);
    cfg.ki = vec1(0.0);
    cfg.kd = vec1(0.0);
    SisoPid pid(cfg);

    auto u = ctrlpp::test::commanded(pid.compute(vec1(100.0), vec1(0.0), Ts));
    REQUIRE_THAT(u[0], WithinAbs(0.0, 1e-15));
}
