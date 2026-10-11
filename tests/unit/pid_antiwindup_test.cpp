#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <cstring>
#include <utility>
#include <algorithm>
#include <type_traits>

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
auto back_calc_refusal(Scalar kp, Scalar ki, Scalar kd, Scalar kb = Scalar{0}) -> ctrlpp::pid_config_error
{
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    using vector_t = typename AwPid::vector_t;
    typename AwPid::config_type cfg{};
    cfg.kp = vector_t::Constant(kp);
    cfg.ki = vector_t::Constant(ki);
    cfg.kd = vector_t::Constant(kd);
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {kb};
    const auto created = AwPid::create(cfg);
    REQUIRE_FALSE(created.has_value());
    return created.error();
}

template <typename Scalar, typename AW, typename Config>
void configure_anti_windup(Config& cfg)
{
    using vector_t = ctrlpp::Vector<Scalar, 1>;
    cfg.kp = vector_t::Constant(Scalar{1});
    cfg.ki = vector_t::Constant(Scalar{1} / Scalar{2});
    if constexpr(std::is_same_v<AW, ctrlpp::anti_windup<ctrlpp::back_calc>>)
        cfg.template policy<AW>().kb = {Scalar{1}};
    if constexpr(std::is_same_v<AW, ctrlpp::anti_windup<ctrlpp::conditional_integration>>)
        cfg.template policy<AW>().error_threshold = {Scalar{2}};
}

struct trajectory_mismatch
{
    std::int32_t cycles;
    std::int32_t first;
};

// Drives a controller with a rate limiter and its twin without one through the same random setpoints,
// which cross zero so that p + (u - p) is inexact on many cycles.
template <typename Scalar, typename AW>
auto rate_limited_mismatch(Scalar rate_max) -> trajectory_mismatch
{
    using Limited = ctrlpp::pid<Scalar, 1, ctrlpp::rate_limit, AW>;
    using Plain = ctrlpp::pid<Scalar, 1, AW>;
    using vector_t = typename Plain::vector_t;
    typename Limited::config_type limited_cfg{};
    typename Plain::config_type plain_cfg{};
    configure_anti_windup<Scalar, AW>(limited_cfg);
    configure_anti_windup<Scalar, AW>(plain_cfg);
    limited_cfg.template policy<ctrlpp::rate_limit>().rate_max = {rate_max};
    auto limited = ctrlpp::test::constructed(Limited::create(limited_cfg));
    auto plain = ctrlpp::test::constructed(Plain::create(plain_cfg));
    std::mt19937 rng{std::mt19937::default_seed};
    std::uniform_real_distribution<Scalar> sp_dist{Scalar{-1}, Scalar{1}};
    std::uniform_real_distribution<Scalar> meas_dist{Scalar{-1} / Scalar{10}, Scalar{1} / Scalar{10}};
    const Scalar dt = Scalar{1} / Scalar{100};
    trajectory_mismatch mismatch{0, -1};
    for(std::int32_t k = 0; k < 20000; ++k)
    {
        const vector_t sp = vector_t::Constant(sp_dist(rng));
        const vector_t meas = vector_t::Constant(meas_dist(rng));
        const auto u_limited = ctrlpp::test::commanded(limited.compute(sp, meas, dt));
        const auto u_plain = ctrlpp::test::commanded(plain.compute(sp, meas, dt));
        if(u_limited[0] == u_plain[0] && limited.integral()[0] == plain.integral()[0])
            continue;
        mismatch.first = mismatch.cycles == 0 ? k : mismatch.first;
        ++mismatch.cycles;
    }
    return mismatch;
}

// The frozen cycle's unconstrained command, kp*1 plus the seeded integrator 1/2, saturates at output_max = 1/2.
template <typename Scalar, typename AW>
auto frozen_through_saturation() -> ctrlpp::pid<Scalar, 1, AW>
{
    using AwPid = ctrlpp::pid<Scalar, 1, AW>;
    using vector_t = typename AwPid::vector_t;
    typename AwPid::config_type cfg{};
    configure_anti_windup<Scalar, AW>(cfg);
    cfg.ki = vector_t::Constant(Scalar{1});
    cfg.output_max = vector_t::Constant(Scalar{1} / Scalar{2});
    auto controller = ctrlpp::test::constructed(AwPid::create(cfg));
    REQUIRE(controller.set_integral(vector_t::Constant(Scalar{1} / Scalar{2})).has_value());
    controller.freeze_integral(true);
    const auto u = ctrlpp::test::commanded(controller.compute(vector_t::Constant(Scalar{1}), vector_t::Zero(), Scalar{1} / Scalar{100}));
    REQUIRE(u[0] == Scalar{1} / Scalar{2});
    return controller;
}

// kp = ki = 1 and dt = 1/64 make every increment e/64, so with dyadic errors each command and integral is exact.
template <typename Controller>
auto unit_gain_config() -> typename Controller::config_type
{
    using vector_t = typename Controller::vector_t;
    using Scalar = typename vector_t::Scalar;
    typename Controller::config_type cfg{};
    cfg.kp = vector_t::Constant(Scalar{1});
    cfg.ki = vector_t::Constant(Scalar{1});
    return cfg;
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
}

TEMPLATE_TEST_CASE("back_calc default kb is refused when the channel has no tracking time",
    "[pid][siso][anti-windup][backcalc][auto-kb]", float, double)
{
    using Scalar = TestType;
    REQUIRE(back_calc_refusal<Scalar>(Scalar{1}, static_cast<Scalar>(-0.1), static_cast<Scalar>(0.1))
        == ctrlpp::pid_config_error::undefined_tracking_time);
    REQUIRE(back_calc_refusal<Scalar>(Scalar{1}, static_cast<Scalar>(-0.5), Scalar{0})
        == ctrlpp::pid_config_error::undefined_tracking_time);
    REQUIRE(back_calc_refusal<Scalar>(Scalar{0}, Scalar{2}, Scalar{0}) == ctrlpp::pid_config_error::undefined_tracking_time);
    REQUIRE(back_calc_refusal<Scalar>(Scalar{0}, Scalar{-2}, Scalar{0}) == ctrlpp::pid_config_error::undefined_tracking_time);

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

    SECTION("a pure-integral channel takes an explicit kb")
    {
        const Scalar ki = 2, meas = -1000;
        const auto expected = reconstructed(Scalar{0}, ki, Scalar{1}, meas);
        REQUIRE(expected.first == Scalar{5});
        REQUIRE(back_calc_first_step<Scalar>(Scalar{0}, ki, Scalar{0}, Scalar{1}, meas) == expected);
    }
}

TEMPLATE_TEST_CASE("back_calc refuses a negative explicit kb", "[pid][siso][anti-windup][backcalc]", float, double)
{
    using Scalar = TestType;
    const Scalar ki = Scalar{1} / Scalar{2};
    for(const Scalar kb : {Scalar{-1}, -std::numeric_limits<Scalar>::denorm_min()})
    {
        CAPTURE(kb);
        REQUIRE(back_calc_refusal<Scalar>(Scalar{1}, ki, Scalar{0}, kb) == ctrlpp::pid_config_error::invalid_back_calculation_gain);
    }
    const Scalar dt = static_cast<Scalar>(0.01), meas = -10;
    const Scalar increment = ki * (Scalar{0} - meas) * dt;
    const Scalar u_raw = Scalar{1} * (Scalar{1} * Scalar{0} - meas) + increment;
    const auto expected = std::pair<Scalar, Scalar>{Scalar{5}, increment + Scalar{1} * (Scalar{5} - u_raw) * dt};
    REQUIRE(back_calc_first_step<Scalar>(Scalar{1}, ki, Scalar{0}, Scalar{1}, meas) == expected);
}

// Without integral action the default kb = sqrt(ki/kd) is sqrt(+-0) = +-0, the limit of an infinite tracking time:
// there is no integrator to drive, so the channel runs exactly as it would without anti-windup.
TEMPLATE_TEST_CASE("back_calc keeps a zero default kb on a mixed-sign channel without integral action",
    "[pid][siso][anti-windup][backcalc][auto-kb]", float, double)
{
    using Scalar = TestType;
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    using PlainPid = ctrlpp::pid<Scalar, 1>;
    using vector_t = typename AwPid::vector_t;
    typename AwPid::config_type cfg{};
    cfg.kp = vector_t::Constant(Scalar{1});
    cfg.kd = vector_t::Constant(-Scalar{1} / Scalar{10});
    cfg.output_min = vector_t::Constant(-Scalar{1} / Scalar{2});
    cfg.output_max = vector_t::Constant(Scalar{1} / Scalar{2});
    typename PlainPid::config_type plain_cfg{};
    plain_cfg.kp = cfg.kp;
    plain_cfg.kd = cfg.kd;
    plain_cfg.output_min = cfg.output_min;
    plain_cfg.output_max = cfg.output_max;
    auto controller = ctrlpp::test::constructed(AwPid::create(cfg));
    auto plain = ctrlpp::test::constructed(PlainPid::create(plain_cfg));
    const Scalar dt = Scalar{1} / Scalar{64};
    for(std::int32_t k = 0; k < 10; ++k)
    {
        CAPTURE(k);
        const auto u = ctrlpp::test::commanded(controller.compute(vector_t::Constant(Scalar{1}), vector_t::Zero(), dt));
        const auto expected = ctrlpp::test::commanded(plain.compute(vector_t::Constant(Scalar{1}), vector_t::Zero(), dt));
        REQUIRE(controller.saturated());
        REQUIRE(std::memcmp(&u[0], &expected[0], sizeof(Scalar)) == 0);
        REQUIRE(std::memcmp(&controller.integral()[0], &plain.integral()[0], sizeof(Scalar)) == 0);
        REQUIRE(controller.integral()[0] == Scalar{0});
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

// rate_max = 1e6 cannot bind: with |sp| <= 1, |meas| <= 1/10 and no saturation, consecutive commands differ by
// at most kp*(2 + 2/10) in the proportional term plus one increment ki*(11/10)*dt, about 2.21, against a
// per-cycle bound rate_max*dt = 1e4.
TEMPLATE_TEST_CASE("an idle rate limiter leaves the command and the integrator bitwise unchanged",
    "[pid][siso][anti-windup][rate-limit]", float, double)
{
    using Scalar = TestType;
    const Scalar idle = static_cast<Scalar>(1e6);

    SECTION("clamping")
    {
        const auto mismatch = rate_limited_mismatch<Scalar, ctrlpp::anti_windup<ctrlpp::clamping>>(idle);
        CAPTURE(mismatch.first);
        REQUIRE(mismatch.cycles == 0);
    }

    SECTION("back_calc")
    {
        const auto mismatch = rate_limited_mismatch<Scalar, ctrlpp::anti_windup<ctrlpp::back_calc>>(idle);
        CAPTURE(mismatch.first);
        REQUIRE(mismatch.cycles == 0);
    }

    SECTION("a binding rate limit changes the trajectory")
    {
        const Scalar binding = Scalar{10};
        REQUIRE(rate_limited_mismatch<Scalar, ctrlpp::anti_windup<ctrlpp::clamping>>(binding).cycles > 0);
        REQUIRE(rate_limited_mismatch<Scalar, ctrlpp::anti_windup<ctrlpp::back_calc>>(binding).cycles > 0);
    }
}

TEMPLATE_TEST_CASE("a frozen integrator holds through saturation under every anti-windup strategy",
    "[pid][siso][anti-windup][freeze]", float, double)
{
    using Scalar = TestType;
    const Scalar half = Scalar{1} / Scalar{2};

    SECTION("back_calc holds while frozen and resumes on unfreeze")
    {
        auto controller = frozen_through_saturation<Scalar, ctrlpp::anti_windup<ctrlpp::back_calc>>();
        REQUIRE(controller.integral()[0] == half);

        controller.freeze_integral(false);
        const Scalar kp = 1, ki = 1, kb = 1, e = 1, dt = Scalar{1} / Scalar{100};
        const auto u = ctrlpp::test::commanded(controller.compute(ctrlpp::Vector<Scalar, 1>::Constant(e), ctrlpp::Vector<Scalar, 1>::Zero(), dt));
        const Scalar integrated = half + ki * e * dt;
        const Scalar u_unconstrained = kp * e + integrated;
        REQUIRE(u[0] == half);
        REQUIRE(controller.integral()[0] == integrated + kb * (half - u_unconstrained) * dt);
    }

    // The frozen cycle had no earlier constraint, so only the freeze keeps its increment out; the next cycle's
    // increment is withheld by the saturation the frozen cycle reported.
    SECTION("clamping holds while frozen and keeps holding against that saturation on unfreeze")
    {
        auto controller = frozen_through_saturation<Scalar, ctrlpp::anti_windup<ctrlpp::clamping>>();
        REQUIRE(controller.integral()[0] == half);

        controller.freeze_integral(false);
        const auto u = ctrlpp::test::commanded(controller.compute(ctrlpp::Vector<Scalar, 1>::Constant(Scalar{1}), ctrlpp::Vector<Scalar, 1>::Zero(), Scalar{1} / Scalar{100}));
        REQUIRE(u[0] == half);
        REQUIRE(controller.integral()[0] == half);
    }

    // |e| = 1 is below the threshold 2, so only the freeze keeps the increment out.
    SECTION("conditional_integration holds while frozen below its threshold")
    {
        REQUIRE(frozen_through_saturation<Scalar, ctrlpp::anti_windup<ctrlpp::conditional_integration>>().integral()[0] == half);
    }
}

// An error of 3/8 adds 3/512 per cycle, so the unconstrained command 3/8 + k*3/512 first passes output_max = 1 on
// cycle 107: 510/512 after cycle 106, 513/512 after cycle 107.
TEMPLATE_TEST_CASE("clamping stores an integrator that reproduces every emitted command",
    "[pid][siso][anti-windup][clamping]", float, double)
{
    using Scalar = TestType;
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;
    using Plain = ctrlpp::pid<Scalar, 1>;
    using vector_t = typename AwPid::vector_t;
    constexpr std::int32_t crossing = 107;
    const Scalar e = Scalar{3} / Scalar{8}, dt = Scalar{1} / Scalar{64}, output_max = 1;
    auto cfg = unit_gain_config<AwPid>();
    auto plain_cfg = unit_gain_config<Plain>();
    cfg.output_max = vector_t::Constant(output_max);
    plain_cfg.output_max = cfg.output_max;
    auto clamped = ctrlpp::test::constructed(AwPid::create(cfg));
    auto plain = ctrlpp::test::constructed(Plain::create(plain_cfg));
    for(std::int32_t k = 1; k <= 200; ++k)
    {
        CAPTURE(k);
        const auto u = ctrlpp::test::commanded(clamped.compute(vector_t::Constant(e), vector_t::Zero(), dt));
        const auto u_plain = ctrlpp::test::commanded(plain.compute(vector_t::Constant(e), vector_t::Zero(), dt));
        REQUIRE(std::min(e + clamped.integral()[0], output_max) == u[0]);
        if(k < crossing)
        {
            REQUIRE(u[0] < output_max);
            REQUIRE(u[0] == u_plain[0]);
            REQUIRE(clamped.integral()[0] == plain.integral()[0]);
        }
        else
            REQUIRE(u[0] == output_max);
    }
}

// Reverse acting: kp = ki = -1 and an error of 3/8 move the command down by 3/512 per cycle onto output_min = -1,
// so the increment, not the error, has the sign of the violated limit.
TEMPLATE_TEST_CASE("clamping holds a reverse-acting integrator once its command saturates",
    "[pid][siso][anti-windup][clamping]", float, double)
{
    using Scalar = TestType;
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;
    using vector_t = typename AwPid::vector_t;
    const Scalar e = Scalar{3} / Scalar{8}, dt = Scalar{1} / Scalar{64}, output_min = -1;
    auto cfg = unit_gain_config<AwPid>();
    cfg.kp = -cfg.kp;
    cfg.ki = -cfg.ki;
    cfg.output_min = vector_t::Constant(output_min);
    auto clamped = ctrlpp::test::constructed(AwPid::create(cfg));
    bool constrained = false;
    Scalar previous = 0;
    for(std::int32_t k = 1; k <= 200; ++k)
    {
        CAPTURE(k);
        const auto u = ctrlpp::test::commanded(clamped.compute(vector_t::Constant(e), vector_t::Zero(), dt));
        const Scalar integral = clamped.integral()[0];
        REQUIRE(std::max(-e + integral, output_min) == u[0]);
        if(constrained)
            REQUIRE(std::abs(integral) <= std::abs(previous));
        constrained = constrained || clamped.saturated();
        previous = integral;
    }
    REQUIRE(constrained);
}

// The error alternates between 1/8, which integrates, and 1/2, past the threshold 1/4; with no output limits the
// command is the proportional term plus the stored integrator on every cycle.
TEMPLATE_TEST_CASE("conditional integration stores an integrator that reproduces every emitted command",
    "[pid][siso][anti-windup][conditional]", float, double)
{
    using Scalar = TestType;
    using AwPid = ctrlpp::pid<Scalar, 1, ctrlpp::anti_windup<ctrlpp::conditional_integration>>;
    using vector_t = typename AwPid::vector_t;
    const Scalar small = Scalar{1} / Scalar{8}, large = Scalar{1} / Scalar{2}, dt = Scalar{1} / Scalar{64};
    auto cfg = unit_gain_config<AwPid>();
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::conditional_integration>>().error_threshold = {Scalar{1} / Scalar{4}};
    auto controller = ctrlpp::test::constructed(AwPid::create(cfg));
    Scalar integrated = 0;
    for(std::int32_t k = 0; k < 64; ++k)
    {
        CAPTURE(k);
        const Scalar e = k % 2 == 0 ? small : large;
        const auto u = ctrlpp::test::commanded(controller.compute(vector_t::Constant(e), vector_t::Zero(), dt));
        integrated += k % 2 == 0 ? small * dt : Scalar{0};
        REQUIRE(u[0] == e + controller.integral()[0]);
        REQUIRE(controller.integral()[0] == integrated);
    }
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

    // Every later cycle follows a cycle saturated at output_max, so each positive increment is withheld.
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
