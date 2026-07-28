// What the oracles in this file decide.
//
// The controller's cycle is fallible and it carries state, so every claim here
// is a claim about the command, the carried state, or the typed rejection that
// protects both:
//
//  * A rejected cycle leaves the integrator and both error histories BITWISE
//    unchanged, and a later valid cycle produces exactly what it would have
//    produced had the bad one never been attempted. Exact, because a rejected
//    cycle performs no arithmetic on the carried state at all.
//  * Where every term of the command is an exact product or an exactly
//    representable sum, the command is asserted with EXACT equality. That covers
//    the proportional-only output, the all-zero-gain output, the extreme-gain
//    outputs, and the known-gain output whose decimal value the arithmetic
//    reproduces bit for bit.
//  * An infinite gain does NOT produce an infinite command: the output
//    saturation clamps against the configured output_max, whose default is the
//    largest finite value of the scalar type. The realized value is that bound,
//    exactly, and asserting it is what proves the clamp runs.
//  * Convergence is asserted against the CLOSED-LOOP MAP of the composite
//    plant-plus-controller system, raised to the number of steps taken and
//    applied to the initial deviation. Neither loop in this file has converged
//    to its setpoint by the step the assertion is made -- one is still 8.2e-4
//    away after ten thousand steps -- so asserting the converged value would
//    manufacture a failure against correct behavior.
//
// What they deliberately do not decide. There is no near-zero-step guard to
// observe: a positive finite step is admitted however small, and the step case
// below asserts that contract rather than the held-output behavior an earlier
// controller had. Nothing here asserts a tuning quality either; the gains are
// given, not designed.

#include "hardening_helpers.h"
#include "ctrlpp/control/pid.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <algorithm>

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
constexpr double eps = std::numeric_limits<double>::epsilon();

// The first-order plant both convergence cases drive.
constexpr double plant_a = 0.9;
constexpr double plant_b = 0.1;

/// Closed-loop map of the plant driven by the position-form controller, in
/// deviation coordinates about the setpoint.
///
/// The composite state is (x_k - r, x_{k-1} - r, I_{k-1} - r), where I is the
/// integrator in output units. Substituting the position form
///     e_k = r - x_k,  I_k = I_{k-1} + ki e_k Ts,
///     d_k = -kd (x_k - x_{k-1}) / Ts,  u_k = kp e_k + I_k + d_k
/// into x_{k+1} = a x_k + b u_k and cancelling the fixed point (r, r, r), which
/// is where e = 0 and the integrator alone supplies u = r, gives the rows below.
/// The derivative term vanishes on the first cycle by construction and also
/// vanishes in this model there, because both measurement histories start equal.
auto closed_loop_map(double kp, double ki, double kd) -> Eigen::Matrix3d
{
    Eigen::Matrix3d M;
    M << plant_a - plant_b * (kp + ki * Ts + kd / Ts), plant_b * (kd / Ts), plant_b,
         1.0, 0.0, 0.0,
         -ki * Ts, 0.0, 1.0;
    return M;
}

auto spectral_radius(Eigen::Matrix3d const& M) -> double
{
    Eigen::EigenSolver<Eigen::Matrix3d> es(M, false);
    double rho = 0.0;
    for(int i = 0; i < 3; ++i)
        rho = std::max(rho, std::abs(es.eigenvalues()(i)));
    return rho;
}

/// Rounded operations along one composite step of plant and controller.
///
/// Enumerated rather than chosen: the error is one product and one difference
/// (two); the proportional term one product (three); the integral term a product
/// by the error, a product by the step, and the accumulation (six); the
/// derivative term a difference, a division by the step, and a product by the
/// gain (nine); the command two additions (eleven); and the plant recursion two
/// products and a sum (fourteen). The reference recursion in the test costs a
/// three-by-three matrix-vector product, three multiplies and two adds per row
/// (fifteen more). Every operation is counted whether or not it rounds, so the
/// count bounds the accumulated error from above rather than describing it.
constexpr int loop_step_ops = 14 + 15;

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

TEST_CASE("PID with an infinite gain saturates at the configured output bound",
          "[pid][hardening][negative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(std::numeric_limits<double>::infinity());
    SisoPid pid(cfg);

    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));

    // An infinite proportional gain against a unit error makes the raw command
    // exactly +infinity, and the output clamp then pins it to output_max, whose
    // default is the largest finite value of the scalar type. That bound is the
    // realized command, exactly. Excluding only NaN passed for every finite
    // value as well, so it could not see either the overflow or the clamp; and
    // the value is NOT infinity, which is what a reader would otherwise assume.
    CHECK(u[0] == std::numeric_limits<double>::max());
    CHECK(u[0] == pid.params().output_max[0]);
}

TEST_CASE("PID admits any positive step and refuses a stopped clock",
          "[pid][hardening][negative]")
{
    // This case previously claimed the controller returns the previous output on
    // a near-zero step. It does not, and it no longer holds a command back at
    // all: a non-positive or non-finite step is REFUSED, and any positive finite
    // step is admitted however small. The old configuration could not have seen
    // either behavior -- it set only the proportional gain, so the command was
    // independent of the step and the assertion held whether or not a guard
    // existed.
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(1.0);
    cfg.kd = vec1(0.5);

    SECTION("a stopped, reversed or garbage clock is refused and mutates nothing")
    {
        SisoPid pid(cfg);
        REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());

        const auto integral_before = pid.integral();
        const auto error_before = pid.error();

        for(double bad_step : {0.0, -Ts, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()})
        {
            CAPTURE(bad_step);
            auto rejected = pid.compute(vec1(1.0), vec1(0.25), bad_step);
            REQUIRE_FALSE(rejected.has_value());
            CHECK(rejected.error() == ctrlpp::pid_step_error::invalid_timestep);
            CHECK(pid.integral() == integral_before);
            CHECK(pid.error() == error_before);
            CHECK(pid.health() == ctrlpp::pid_health::ok);
        }
    }

    SECTION("a step below the resolution of the clock is a valid step")
    {
        constexpr double tiny_step = 1e-15;
        SisoPid pid(cfg);

        auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
        auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.25), tiny_step));

        // The command is the full position form at that step, and the derivative
        // term dominates it because the measurement moved while the clock barely
        // did: the difference quotient is divided by the step. That is the
        // property a held output would destroy, and it is what makes the step
        // observable at all.
        //
        // Step one: e = 1, so the proportional term is 2 and the integrator takes
        // ki e Ts. The derivative term is zero on the first cycle by
        // construction. Step two: e = 0.75, the integrator takes a further
        // ki e dt, and the derivative term is -kd (0.25 - 0) / dt.
        double const integral_1 = 1.0 * 1.0 * Ts;
        double const integral_2 = integral_1 + 1.0 * 0.75 * tiny_step;
        double const derivative_2 = -0.5 * (0.25 - 0.0) / tiny_step;
        double const expected_2 = 2.0 * 0.75 + integral_2 + derivative_2;

        // Four roundings stand behind the derivative term and the sum that
        // carries it: the difference, the division, the gain product and the
        // final accumulation. The scale is the term that entered largest.
        constexpr int step_command_ops = 4;
        CAPTURE(u1[0], u2[0], expected_2);
        REQUIRE(std::abs(u2[0] - expected_2)
                <= step_command_ops * eps * std::abs(derivative_2));

        // And it is emphatically not the previous command.
        REQUIRE(std::abs(u2[0] - u1[0]) > std::abs(derivative_2) / 2.0);
    }
}

TEST_CASE("PID P-only exact output", "[pid][hardening][precision]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    SisoPid pid(cfg);

    // A unit gain against a unit error, with both other terms an exact product of
    // zero: the command is exactly one and the case name already said so.
    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE(u[0] == 1.0);
}

TEST_CASE("PID precision with known PID gains", "[pid][hardening][precision]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    cfg.kd = vec1(0.0);
    SisoPid pid(cfg);

    // Step 1: e = 1, P = 2, I = 0.5 * 1 * 0.01, D = 0.
    //
    // Exact, and the reason is worth stating because 0.01 is not representable.
    // The integral increment is that inexact step halved, and halving is exact in
    // a binary radix, so the increment is the nearest double to 0.005 -- which is
    // also what the decimal literal 0.005 denotes, by the same scaling argument.
    // The remaining sum 2 + 0.005 rounds to the nearest double to 2.005. Every
    // step of the chain lands on the literal below, bit for bit.
    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE(u[0] == 2.005);
}

namespace {

/// Drive the plant with the controller for `steps` cycles from rest, and assert
/// the realized state against the closed-loop map raised to the same power.
///
/// The composite system is linear, so the deviation from the setpoint after
/// `steps` cycles is exactly the map applied `steps` times to the initial
/// deviation. That is the property "stable gains" names, and it fails an
/// implementation that diverges, that settles on the wrong value, or that
/// converges at a rate its own gains do not produce -- none of which a one-sided
/// magnitude threshold can see.
///
/// The rounding budget is geometric, not additive: a perturbation injected at
/// cycle k is itself contracted by the remaining cycles, so the accumulated
/// deviation is bounded by the per-step count times the machine epsilon times the
/// scale that entered, summed over a geometric series with ratio equal to the
/// closed-loop spectral radius.
void assert_closed_loop_trajectory(double kp, double ki, double kd, int steps)
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(kp);
    cfg.ki = vec1(ki);
    cfg.kd = vec1(kd);
    SisoPid pid(cfg);

    constexpr double setpoint = 1.0;
    double x = 0.0;
    double scale = setpoint;

    for(int k = 0; k < steps; ++k)
    {
        auto u = ctrlpp::test::commanded(pid.compute(vec1(setpoint), vec1(x), Ts));
        scale = std::max({scale, std::abs(u[0]), std::abs(x)});
        x = plant_a * x + plant_b * u[0];
    }

    auto const M = closed_loop_map(kp, ki, kd);
    double const rho = spectral_radius(M);
    REQUIRE(rho < 1.0);

    Eigen::Vector3d z(-setpoint, -setpoint, -setpoint);
    for(int k = 0; k < steps; ++k)
        z = (M * z).eval();

    double const budget = loop_step_ops * eps * scale / (1.0 - rho);

    CAPTURE(steps, rho, x, z(0), scale, budget);
    REQUIRE(std::abs((x - setpoint) - z(0)) <= budget);
}

}

TEST_CASE("PID stays bounded over 10000 steps with stable gains", "[pid][hardening][stability]")
{
    // Measured, and recorded because the intuition is wrong: with these gains the
    // loop has NOT converged after ten thousand cycles. The integral gain is
    // small enough that the remaining deviation is still 8.2e-4, so asserting the
    // setpoint would fail against entirely correct behavior. What the closed-loop
    // map asserts instead is the whole transient, including the part still in
    // flight.
    assert_closed_loop_trajectory(0.5, 0.1, 0.01, 10000);
}

TEST_CASE("PID with integral action converges to zero error", "[pid][hardening][convergence]")
{
    // Same form with a decade more integral action: the deviation after five
    // thousand cycles is 1.6e-8, which is still eight decades above the machine
    // floor and is exactly what the map predicts.
    assert_closed_loop_trajectory(2.0, 1.0, 0.0, 5000);
}

TEST_CASE("PID extreme gains produce exactly the commanded sum", "[pid][hardening][robustness]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1e12);
    cfg.ki = vec1(1e12);
    SisoPid pid(cfg);

    // Both terms are exact: 1e12 and 1e12 * 0.01 are representable, and so is
    // their sum. Finiteness held for every value the controller could return,
    // including a wrong one, whereas these are the only values it may return.
    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    CHECK(u[0] == 1e12 + 1e12 * Ts);

    // Second cycle: the error halves, the proportional term with it, and the
    // integrator accumulates the second increment on top of the first.
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.5), Ts));
    CHECK(u2[0] == 1e12 * 0.5 + (1e12 * Ts + 1e12 * 0.5 * Ts));
}

TEST_CASE("PID zero gains produce zero output", "[pid][hardening][robustness]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(0.0);
    cfg.ki = vec1(0.0);
    cfg.kd = vec1(0.0);
    SisoPid pid(cfg);

    // Every term is an exact product of a zero gain, so the sum is exactly zero
    // and a tolerance said less than the arithmetic guarantees.
    auto u = ctrlpp::test::commanded(pid.compute(vec1(100.0), vec1(0.0), Ts));
    REQUIRE(u[0] == 0.0);
}
