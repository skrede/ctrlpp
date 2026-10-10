#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

// Bound for comparing a value the test reconstructs against the same value the
// controller computes with a different operation order. Each of the arithmetic
// operations composing the stored integral is rounded once in the controller and
// once in the test, so the accumulated forward error is bounded by twice the
// operation count times one unit in the last place at the working scale.
double integral_reconstruction_tol(double scale)
{
    constexpr int arithmetic_ops = 10; // products/sums forming the stored integral
    return 2.0 * arithmetic_ops * std::numeric_limits<double>::epsilon() * std::abs(scale);
}

template <typename Scalar>
auto back_calc_first_step(Scalar kp, Scalar ki, Scalar kd, Scalar kb_config, Scalar meas) -> std::pair<Scalar, Scalar>
{
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    using vector_t = typename AwPid::vector_t;
    typename AwPid::config_type cfg{};
    cfg.kp = vector_t::Constant(kp);
    cfg.ki = vector_t::Constant(ki);
    cfg.kd = vector_t::Constant(kd);
    cfg.output_max = vector_t::Constant(Scalar{5});
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {kb_config};
    auto controller = ctrlpp::test::constructed(AwPid::create(cfg));
    const Scalar dt = static_cast<Scalar>(0.01);
    const auto u = ctrlpp::test::commanded(controller.compute(vector_t::Zero(), vector_t::Constant(meas), dt));
    return {u[0], controller.integral()[0]};
}

template <typename Scalar>
auto back_calc_refusal(Scalar kp, Scalar ki, Scalar kd) -> ctrlpp::pid_config_error
{
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    using vector_t = typename AwPid::vector_t;
    typename AwPid::config_type cfg{};
    cfg.kp = vector_t::Constant(kp);
    cfg.ki = vector_t::Constant(ki);
    cfg.kd = vector_t::Constant(kd);
    const auto created = AwPid::create(cfg);
    REQUIRE_FALSE(created.has_value());
    return created.error();
}

}

TEST_CASE("Without anti_windup: integral winds up unboundedly", "[pid][siso][windup]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    auto pid = ctrlpp::test::constructed(SisoPid::create(cfg));

    // Constant error of 1.0 for many steps
    for (int i = 0; i < 1000; ++i)
        REQUIRE(pid.compute(vec1(10.0), vec1(0.0), Ts).has_value());

    // Integral should be large (no anti-windup to stop it)
    // I = Ki * e * Ts * 1000 = 1 * 10 * 0.01 * 1000 ~ 100
    REQUIRE_THAT(pid.integral()[0], WithinAbs(100.0, 1e-6));
}

TEST_CASE("back_calc anti-windup limits integral growth during saturation",
    "[pid][siso][anti-windup][backcalc]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::back_calc>;
    using AwPid = ctrlpp::pid<double, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    cfg.template policy<AW>().kb = {1.0};
    auto pid = ctrlpp::test::constructed(AwPid::create(cfg));

    // Run many steps with large error to cause saturation
    for (int i = 0; i < 1000; ++i)
        REQUIRE(pid.compute(vec1(10.0), vec1(0.0), Ts).has_value());

    // Integral should be bounded (back-calc feedback limits growth)
    // Without anti-windup it would be 100.0
    REQUIRE(pid.integral()[0] < 50.0);
}

// The back-calculation gain kb multiplies the saturation error (u_sat - u_raw, in
// output units) into the integrator, whose rate is in output-per-time, so kb must
// carry units of 1/time. The dimensionally-correct auto default is the Astrom
// tracking-time-constant form kb = 1/Tt with Tt = sqrt(Ti*Td), which in the internal
// parallel gains is kb = sqrt(ki/kd); with no derivative action it collapses to the
// integral-time reciprocal kb = ki/kp. Each case below drives one step from reset and
// reconstructs the stored integral from the model, which pins the exact kb the
// controller used (the first step has no derivative term, so the integral after the
// step is ki*e*dt from integration plus kb*(u_sat - u_raw)*dt from back-calculation).
TEST_CASE("back_calc default Kb auto-computation", "[pid][siso][anti-windup][backcalc][auto-kb]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::back_calc>;

    const double sp = 10.0;
    const double meas = 0.0;
    const double e = sp - meas; // setpoint weight b defaults to 1

    SECTION("kb = sqrt(ki/kd) when kd != 0") {
        using AwPid = ctrlpp::pid<double, 1, AW>;
        const double kp = 2.0, ki = 8.0, kd = 2.0, out_max = 5.0;

        AwPid::config_type cfg{};
        cfg.kp = vec1(kp);
        cfg.ki = vec1(ki);
        cfg.kd = vec1(kd);
        cfg.output_max = vec1(out_max);
        // kb left at 0 -> auto default. sqrt(ki/kd) = sqrt(8/2) = 2, which is
        // distinct from the dimensionally-wrong sqrt(ki*kd) = sqrt(16) = 4.
        auto pid = ctrlpp::test::constructed(AwPid::create(cfg));

        auto u = ctrlpp::test::commanded(pid.compute(vec1(sp), vec1(meas), Ts));

        const double increment = ki * e * Ts;
        const double u_raw = kp * e + increment; // first step: no derivative, no ff
        const double u_sat = std::min(u_raw, out_max);
        const double kb = std::sqrt(ki / kd);
        const double integral_expected = increment + kb * (u_sat - u_raw) * Ts;

        REQUIRE_THAT(u[0], WithinAbs(u_sat, integral_reconstruction_tol(u_raw)));
        REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_expected, integral_reconstruction_tol(u_raw)));
    }

    SECTION("kb = ki/kp fallback when kd == 0") {
        using AwPid = ctrlpp::pid<double, 1, AW>;
        const double kp = 4.0, ki = 2.0, out_max = 5.0;

        AwPid::config_type cfg{};
        cfg.kp = vec1(kp);
        cfg.ki = vec1(ki);
        // kd left at 0 -> auto default kb = ki/kp = 0.5, distinct from the old kb = ki = 2.
        cfg.output_max = vec1(out_max);
        auto pid = ctrlpp::test::constructed(AwPid::create(cfg));

        auto u = ctrlpp::test::commanded(pid.compute(vec1(sp), vec1(meas), Ts));

        const double increment = ki * e * Ts;
        const double u_raw = kp * e + increment;
        const double u_sat = std::min(u_raw, out_max);
        const double kb = ki / kp;
        const double integral_expected = increment + kb * (u_sat - u_raw) * Ts;

        REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_expected, integral_reconstruction_tol(u_raw)));
    }

    SECTION("pure-I controller disables back-calc (kb = 0) with kp == 0 and kd == 0") {
        using AwPid = ctrlpp::pid<double, 1, AW>;
        const double ki = 2.0, out_max = 0.1; // kp = kd = 0

        AwPid::config_type cfg{};
        cfg.ki = vec1(ki);
        cfg.output_max = vec1(out_max);
        auto pid = ctrlpp::test::constructed(AwPid::create(cfg));

        auto u = ctrlpp::test::commanded(pid.compute(vec1(sp), vec1(meas), Ts));

        // No proportional or derivative reference for a tracking time, so kb defaults
        // to 0: the output still saturates, but the integrator receives no
        // back-calculation feedback and holds only its integration increment.
        const double increment = ki * e * Ts;
        REQUIRE(u[0] == out_max); // saturated
        REQUIRE_THAT(pid.integral()[0], WithinAbs(increment, integral_reconstruction_tol(increment)));
    }
}

TEMPLATE_TEST_CASE("back_calc default kb is refused when the gains do not share a sign",
    "[pid][siso][anti-windup][backcalc][auto-kb]", float, double)
{
    using Scalar = TestType;
    REQUIRE(back_calc_refusal<Scalar>(Scalar{1}, static_cast<Scalar>(-0.1), static_cast<Scalar>(0.1))
        == ctrlpp::pid_config_error::undefined_tracking_time);
    REQUIRE(back_calc_refusal<Scalar>(Scalar{1}, static_cast<Scalar>(-0.5), Scalar{0})
        == ctrlpp::pid_config_error::undefined_tracking_time);

    // The first step saturates at output_max = 5 from setpoint 0; the expected command and integral are
    // formed in the controller's operation order, so they match bit for bit.
    const Scalar dt = static_cast<Scalar>(0.01);
    const Scalar sp = 0;
    const auto reconstructed = [&](Scalar kp, Scalar ki, Scalar kb, Scalar meas) {
        const Scalar increment = ki * (sp - meas) * dt;
        const Scalar u_raw = kp * (Scalar{1} * sp - meas) + increment;
        const Scalar u_sat = std::min(u_raw, Scalar{5});
        return std::pair<Scalar, Scalar>{u_sat, increment + kb * (u_sat - u_raw) * dt};
    };

    SECTION("reverse-acting gains of one sign take kb = sqrt(ki/kd)")
    {
        const Scalar kp = -1, ki = static_cast<Scalar>(-0.5), kd = static_cast<Scalar>(-0.1), meas = 10;
        const auto expected = reconstructed(kp, ki, std::sqrt(ki / kd), meas);
        REQUIRE(expected.first == Scalar{5});
        REQUIRE(back_calc_first_step<Scalar>(kp, ki, kd, Scalar{0}, meas) == expected);
    }

    SECTION("an explicit kb is accepted whatever the gain signs")
    {
        const Scalar kp = 1, ki = static_cast<Scalar>(-0.1), kd = static_cast<Scalar>(0.1), meas = -10;
        const auto expected = reconstructed(kp, ki, Scalar{2}, meas);
        REQUIRE(expected.first == Scalar{5});
        REQUIRE(back_calc_first_step<Scalar>(kp, ki, kd, Scalar{2}, meas) == expected);
    }
}

TEST_CASE("back_calc anti-windup limits windup during a rate-limited ramp",
    "[pid][siso][anti-windup][backcalc][rate-limit]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::back_calc>;
    using RlPid = ctrlpp::pid<double, 1, AW, ctrlpp::rate_limit>;

    RlPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.template policy<AW>().kb = {1.0};
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {1.0}; // slow actuator, no output clamp

    auto pid = ctrlpp::test::constructed(RlPid::create(cfg));

    // Large constant setpoint: the unconstrained command wants to jump far, but the
    // rate limiter allows only rate_max*dt per step. Because back-calculation feeds back
    // the difference against the unconstrained command (not the already rate-limited
    // value), it sees the rate-limit constraint and holds the integrator.
    for (int i = 0; i < 1000; ++i)
        REQUIRE(pid.compute(vec1(10.0), vec1(0.0), Ts).has_value());

    // Feeding back only the rate-limited value would make the saturation error zero, so
    // the integral would wind up to ki*e*dt*steps = 100. The rate-aware feedback keeps it
    // bounded.
    REQUIRE(pid.integral()[0] < 50.0);
}

TEST_CASE("clamping anti-windup freezes integral during saturation",
    "[pid][siso][anti-windup][clamping]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::clamping>;
    using AwPid = ctrlpp::pid<double, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.output_max = vec1(5.0);
    auto pid = ctrlpp::test::constructed(AwPid::create(cfg));

    // Step until saturation
    double integral_at_saturation = 0.0;
    bool found_saturation = false;
    for (int i = 0; i < 100; ++i) {
        REQUIRE(pid.compute(vec1(10.0), vec1(0.0), Ts).has_value());
        if (pid.saturated() && !found_saturation) {
            found_saturation = true;
            integral_at_saturation = pid.integral()[0];
        }
    }
    REQUIRE(found_saturation);

    // After many more steps during saturation, integral should be frozen
    // (error > 0 and integral > 0 during saturation -> undo increment)
    for (int i = 0; i < 100; ++i)
        REQUIRE(pid.compute(vec1(10.0), vec1(0.0), Ts).has_value());

    // Integral should have stayed near the saturation point
    REQUIRE_THAT(pid.integral()[0], WithinAbs(integral_at_saturation, tol));
}

TEST_CASE("conditional_integration freezes integral when error exceeds threshold",
    "[pid][siso][anti-windup][conditional]")
{
    using AW = ctrlpp::anti_windup<ctrlpp::conditional_integration>;
    using AwPid = ctrlpp::pid<double, 1, AW>;

    AwPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(1.0);
    cfg.template policy<AW>().error_threshold = {2.0};
    auto pid = ctrlpp::test::constructed(AwPid::create(cfg));

    // Error = 5 > threshold 2 -> integral should not accumulate
    for (int i = 0; i < 100; ++i)
        REQUIRE(pid.compute(vec1(5.0), vec1(0.0), Ts).has_value());

    REQUIRE_THAT(pid.integral()[0], WithinAbs(0.0, tol));

    // Error = 1 < threshold 2 -> integral should accumulate
    pid.reset();
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());
    REQUIRE_THAT(pid.integral()[0], WithinAbs(1.0 * 1.0 * Ts, tol));
}
