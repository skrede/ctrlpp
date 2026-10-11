#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

template <typename S>
auto vs(S v) -> ctrlpp::Vector<S, 1>
{
    return ctrlpp::Vector<S, 1>::Constant(v);
}

// Inputs lie on a 2^-3 grid below 2^2 and the gains and steps are powers of two, so over eight steps
// every intermediate of either form, including the halving of the alpha = 1/2 filters, lies on a
// 2^-18 grid below 2^4: 22 bits, inside float's 24-bit significand. Both sides are then exact and
// the running sum equals the position output bit for bit; no tolerance is needed.
template <typename S, typename... Policies>
auto velocity_position_mismatches(const auto& configure_velocity, const auto& configure_position, bool alternate_dt) -> std::int32_t
{
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form, Policies...>;
    using PPid = ctrlpp::pid<S, 1, Policies...>;
    constexpr std::array<double, 8> setpoints{1.0, 1.0, 1.5, 1.5, 2.0, 0.5, 0.5, 1.0};
    constexpr std::array<double, 8> measurements{0.0, 0.125, 0.375, 0.625, 1.0, 1.125, 0.875, 0.75};
    const auto gains = [](auto& cfg) { cfg.kp = vs(S{1}); cfg.ki = vs(S{1} / S{2}); cfg.kd = vs(S{1} / S{8}); };
    typename VPid::config_type vcfg{};
    typename PPid::config_type pcfg{};
    gains(vcfg);
    gains(pcfg);
    configure_velocity(vcfg);
    configure_position(pcfg);
    auto velocity = ctrlpp::test::constructed(VPid::create(vcfg));
    auto position = ctrlpp::test::constructed(PPid::create(pcfg));
    S sum{0};
    std::int32_t mismatches = 0;
    for(std::size_t k = 0; k < setpoints.size(); ++k)
    {
        const S dt = alternate_dt && k % 2 == 1 ? S{1} / S{32} : S{1} / S{64};
        const auto sp = vs(static_cast<S>(setpoints[k]));
        const auto meas = vs(static_cast<S>(measurements[k]));
        sum += ctrlpp::test::commanded(velocity.compute(sp, meas, dt))[0];
        mismatches += sum == ctrlpp::test::commanded(position.compute(sp, meas, dt))[0] ? 0 : 1;
    }
    return mismatches;
}

}

TEST_CASE("velocity_form P-only: delta_u = Kp*(e(k)-e(k-1))", "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    // Step 1: e=1.0, prev_e=0.0
    // delta_u = Kp*(1-0) + 0 + 0 = 2.0
    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // Step 2: e=1.0, prev_e=1.0
    // delta_u = Kp*(1-1) = 0.0
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u2[0], WithinAbs(0.0, tol));

    // Step 3: e=0.5, prev_e=1.0
    // delta_u = Kp*(0.5-1.0) = -1.0
    auto u3 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.5), Ts));
    REQUIRE_THAT(u3[0], WithinAbs(-1.0, tol));
}

TEST_CASE("velocity_form PI: includes Ki*e*dt incremental term", "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.5);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    // Step 1: e=1.0, prev_e=0
    // dP = 1*(1-0) = 1.0, dI = 0.5*1*0.01 = 0.005
    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    // Step 2: e=1.0, prev_e=1.0
    // dP = 1*(1-1) = 0, dI = 0.5*1*0.01 = 0.005
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u2[0], WithinAbs(0.005, tol));

    // Step 3: e=0.5
    // dP = 1*(0.5-1) = -0.5, dI = 0.5*0.5*0.01 = 0.0025
    auto u3 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.5), Ts));
    REQUIRE_THAT(u3[0], WithinAbs(-0.5 + 0.0025, tol));
}

TEMPLATE_TEST_CASE("velocity_form PID emits no derivative increment on the first step, then a first and a second difference",
    "[pid][siso][velocity-form]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    const S kp{2};
    const S ki = S{1} / S{2};
    const S kd = S{1} / S{8};
    const S dt = S{1} / S{64};
    typename VPid::config_type cfg{};
    cfg.kp = vs(kp);
    cfg.ki = vs(ki);
    cfg.kd = vs(kd);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    const S p1 = kp * (S{1} - S{0});
    const S du1 = ((p1 - S{0}) + (ki * S{1}) * dt) + (S{0} - S{0});
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == du1);

    const S m2 = S{1} / S{4};
    const S p2 = kp * (S{1} - m2);
    const S d2 = -(kd * ((m2 - S{0}) / dt));
    const S du2 = ((p2 - p1) + (ki * (S{1} - m2)) * dt) + (d2 - S{0});
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(m2), dt))[0] == du2);

    const S d3 = -(kd * ((m2 - m2) / dt));
    const S du3 = ((p2 - p2) + (ki * (S{1} - m2)) * dt) + (d3 - d2);
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(m2), dt))[0] == du3);
}

TEMPLATE_TEST_CASE("velocity_form PD setpoint step emits kp*e once, with no derivative kick and no parked offset",
    "[pid][siso][velocity-form]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    const S dt = S{1} / S{1000};
    for(const bool on_error : {false, true})
    {
        CAPTURE(on_error);
        typename VPid::config_type cfg{};
        cfg.kp = vs(S{1});
        cfg.kd = vs(S{1});
        cfg.output_min = vs(S{-100});
        cfg.output_max = vs(S{100});
        cfg.derivative_on_error = on_error;
        auto pid = ctrlpp::test::constructed(VPid::create(cfg));
        for(std::int32_t step = 1; step <= 10; ++step)
        {
            CAPTURE(step);
            CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == (step == 1 ? S{1} : S{0}));
            CHECK_FALSE(pid.saturated());
        }
    }
}

TEST_CASE("velocity_form steady-state: delta_u converges to Ki*e*dt",
    "[pid][siso][velocity-form]")
{
    using VPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(2.0);
    cfg.kd = vec1(0.1);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    // First few steps build up history, then with constant error:
    // dP = 0 (error unchanged), dD = 0 (measurement unchanged), dI = Ki*e*dt
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());  // step 1
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());  // step 2
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());  // step 3

    // From step 3 onward, e is constant at 1.0
    // delta_u should be Ki*e*dt = 2*1*0.01 = 0.02
    auto u4 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u4[0], WithinAbs(0.02, tol));

    auto u5 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u5[0], WithinAbs(0.02, tol));
}

TEST_CASE("velocity_form reset clears history", "[pid][siso][velocity-form][reset]")
{
    using VPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.ki = vec1(0.5);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    // Build up history
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), Ts).has_value());

    pid.reset();

    // After reset: behaves as the first step
    auto u = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    // dP = 2*(1-0) = 2, dI = 0.5*1*0.01 = 0.005, dD = 0 (Kd=0)
    REQUIRE_THAT(u[0], WithinAbs(2.005, tol));
}

TEST_CASE("velocity_form injects the change in feed-forward so a constant feed-forward does not drift",
    "[pid][siso][velocity-form][feed-forward]")
{
    auto ff_func = [](const Vec1& sp, double) -> Vec1 { return sp * 0.5; };
    using FF = ctrlpp::feed_forward<decltype(ff_func)>;
    using VFfPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form, FF>;

    VFfPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.ki = vec1(0.0);
    cfg.template policy<FF>().ff_func = ff_func;
    auto pid = ctrlpp::test::constructed(VFfPid::create(cfg));

    // Step 1: dP = 1*(1-0) = 1.0. The feed-forward level rises from 0 to sp*0.5 = 0.5,
    // so its increment is 0.5. delta_u = 1.0 + 0.5 = 1.5.
    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(1.5, tol));

    // Step 2 onward: dP = 0 and the feed-forward level is unchanged (constant sp), so
    // the feed-forward increment is 0 and it contributes nothing to delta_u. Injecting
    // the absolute feed-forward instead would emit 0.5 every step and the actuator would
    // drift without bound.
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u2[0], WithinAbs(0.0, tol));

    auto u3 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u3[0], WithinAbs(0.0, tol));
}

TEMPLATE_TEST_CASE("velocity_form tracking re-seeds the accumulated output to the clamped tracking signal",
    "[pid][siso][velocity-form][tracking]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    typename VPid::config_type cfg{};
    cfg.kp = vs(S{1});
    cfg.output_min = vs(S{-10});
    cfg.output_max = vs(S{10});
    const S dt = S{1} / S{64};
    const auto tracked_then_stepped = [&](S w) {
        auto pid = ctrlpp::test::constructed(VPid::create(cfg));
        CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{1});
        CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt, vs(w)))[0] == S{0});
        CHECK(pid.integral()[0] == S{0});
        const S increment = ctrlpp::test::commanded(pid.compute(vs(S{3}), vs(S{0}), dt))[0];
        return std::pair{increment, pid.saturated()};
    };

    SECTION("a target inside the limits")
    {
        const auto [increment, saturated] = tracked_then_stepped(S{19} / S{2});
        CHECK(increment == S{1} / S{2});
        CHECK(saturated);
    }
    SECTION("a target past the upper limit re-seeds to the limit")
    {
        const auto [increment, saturated] = tracked_then_stepped(S{50});
        CHECK(increment == S{0});
        CHECK(saturated);
    }
}

TEST_CASE("velocity_form clamps the accumulated output so it can rise and fall under an asymmetric limit",
    "[pid][siso][velocity-form][clamping]")
{
    using VPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(10.0);
    cfg.output_min = vec1(0.0);   // valve-like: the output cannot go below zero
    cfg.output_max = vec1(1.0);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    // The output limits bound the accumulated output, not the per-step increment, so
    // the running sum of the emitted increments must stay in [output_min, output_max].
    double accumulated = 0.0;

    // Step 1: e=1, dP=10 -> target accumulated 0+10=10, clamped to output_max=1,
    // so the emitted increment is 1.0.
    double d1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts))[0];
    accumulated += d1;
    REQUIRE_THAT(accumulated, WithinAbs(1.0, tol));

    // Step 2: e=1, prev_e=1, dP=0 -> accumulated stays at output_max, increment 0.
    double d2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts))[0];
    accumulated += d2;
    REQUIRE_THAT(accumulated, WithinAbs(1.0, tol));

    // Step 3: error flips negative, dP=10*(-1-1)=-20 -> target 1-20=-19, clamped to
    // output_min=0. The emitted increment is NEGATIVE (-1.0): the output falls back to
    // its lower limit. The old increment clamp (delta_u >= output_min = 0) forbade this,
    // pinning the output at output_max.
    double d3 = ctrlpp::test::commanded(pid.compute(vec1(0.0), vec1(1.0), Ts))[0];
    accumulated += d3;
    REQUIRE(d3 < 0.0);
    REQUIRE_THAT(accumulated, WithinAbs(0.0, tol));
}

TEST_CASE("velocity_form zero dt is rejected rather than answered with the held output",
    "[pid][siso][velocity-form][edge-case]")
{
    using VPid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    VPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(2.0, tol));

    // A stopped clock names the fault instead of handing back the stored
    // output, which the caller could not tell apart from a computed one. The
    // accumulated output and both error histories are left untouched.
    const auto error_before = pid.error();
    auto stopped = pid.compute(vec1(5.0), vec1(0.0), 0.0);
    REQUIRE_FALSE(stopped.has_value());
    CHECK(stopped.error() == ctrlpp::pid_step_error::invalid_timestep);
    CHECK(pid.error() == error_before);
    CHECK(pid.health() == ctrlpp::pid_health::ok);
}

TEMPLATE_TEST_CASE("velocity_form reports saturated exactly on cycles whose accumulated output is clamped",
    "[pid][siso][velocity-form][clamping]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    typename VPid::config_type cfg{};
    cfg.kp = vs(S{10});
    cfg.output_min = vs(S{0});
    cfg.output_max = vs(S{1});
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));
    const S dt = S{1} / S{64};

    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{1});
    CHECK(pid.saturated());
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{0});
    CHECK_FALSE(pid.saturated());
}

TEMPLATE_TEST_CASE("velocity_form derivative increment is exact under non-uniform sampling",
    "[pid][siso][velocity-form][derivative]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    const S kd = S{1} / S{8};
    for(const bool uniform : {false, true})
    {
        CAPTURE(uniform);
        typename VPid::config_type cfg{};
        cfg.kd = vs(kd);
        auto pid = ctrlpp::test::constructed(VPid::create(cfg));
        S t{0};
        for(std::int32_t step = 1; step <= 8; ++step)
        {
            CAPTURE(step);
            const S dt = uniform || step % 2 == 1 ? S{1} / S{64} : S{1} / S{32};
            t += dt;
            const auto du = ctrlpp::test::commanded(pid.compute(vs(S{0}), vs(t), dt));
            CHECK(du[0] == (step == 2 ? -kd : S{0}));
        }
    }
}

TEMPLATE_TEST_CASE("velocity_form freeze_integral zeroes the integral increment until unfrozen",
    "[pid][siso][velocity-form][freeze]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    const S ki = S{1} / S{2};
    const S dt = S{1} / S{64};
    typename VPid::config_type cfg{};
    cfg.kp = vs(S{1});
    cfg.ki = vs(ki);
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));

    pid.freeze_integral(true);
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{1});
    for(std::int32_t step = 2; step <= 4; ++step)
        CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{0});
    pid.freeze_integral(false);
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == (ki * S{1}) * dt);
}

TEMPLATE_TEST_CASE("velocity_form running sum equals the position form output when no limit binds",
    "[pid][siso][velocity-form][differential]", float, double)
{
    using S = TestType;
    const auto base = [](auto&) {};
    SECTION("setpoint weights with the derivative on the weighted error")
    {
        const auto weights = [](auto& cfg) { cfg.b = vs(S{1} / S{2}); cfg.c = vs(S{1} / S{4}); cfg.derivative_on_error = true; };
        CHECK(velocity_position_mismatches<S>(weights, weights, false) == 0);
    }
    SECTION("derivative filter")
    {
        const auto filter = [](auto& cfg) { cfg.template policy<ctrlpp::deriv_filter>().n = {S{8}}; };
        CHECK(velocity_position_mismatches<S, ctrlpp::deriv_filter>(filter, filter, false) == 0);
    }
    SECTION("forward Euler integration")
    {
        CHECK(velocity_position_mismatches<S, ctrlpp::forward_euler>(base, base, false) == 0);
    }
    SECTION("Tustin integration")
    {
        CHECK(velocity_position_mismatches<S, ctrlpp::tustin>(base, base, false) == 0);
    }
    SECTION("setpoint and measurement filters")
    {
        const auto filters = [](auto& cfg) {
            cfg.template policy<ctrlpp::setpoint_filter>().tf = {S{1} / S{64}};
            cfg.template policy<ctrlpp::pv_filter>().tf = {S{1} / S{64}};
        };
        CHECK(velocity_position_mismatches<S, ctrlpp::setpoint_filter, ctrlpp::pv_filter>(filters, filters, false) == 0);
    }
    SECTION("feed-forward")
    {
        auto ff_func = [](const ctrlpp::Vector<S, 1>& sp, S) -> ctrlpp::Vector<S, 1> { return sp / S{4}; };
        using FF = ctrlpp::feed_forward<decltype(ff_func)>;
        const auto feed = [&](auto& cfg) { cfg.template policy<FF>().ff_func = ff_func; };
        CHECK(velocity_position_mismatches<S, FF>(feed, feed, false) == 0);
    }
    SECTION("non-uniform sampling")
    {
        CHECK(velocity_position_mismatches<S>(base, base, true) == 0);
    }
    SECTION("negative control: a different setpoint weight on one side is detected")
    {
        const auto half = [](auto& cfg) { cfg.b = vs(S{1} / S{2}); cfg.c = vs(S{1} / S{4}); cfg.derivative_on_error = true; };
        const auto unit = [](auto& cfg) { cfg.c = vs(S{1} / S{4}); cfg.derivative_on_error = true; };
        CHECK(velocity_position_mismatches<S>(half, unit, false) > 0);
    }
}

TEMPLATE_TEST_CASE("velocity_form clips each increment to the rate bound and carries no remainder",
    "[pid][siso][velocity-form][rate-limit]", float, double)
{
    using S = TestType;
    using VPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form, ctrlpp::rate_limit>;
    typename VPid::config_type cfg{};
    cfg.kp = vs(S{10});
    cfg.output_min = vs(S{-100});
    cfg.output_max = vs(S{100});
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {S{1}};
    auto pid = ctrlpp::test::constructed(VPid::create(cfg));
    const S dt = S{1} / S{64};

    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{1} / S{64});
    CHECK_FALSE(pid.saturated());
    CHECK(ctrlpp::test::commanded(pid.compute(vs(S{1}), vs(S{0}), dt))[0] == S{0});
    CHECK_FALSE(pid.saturated());
}

// rate_max = 1e6 cannot bind: |dP| <= 2 for these inputs, |D| <= kd*(3/8)/dt = 3 so |dD| <= 6, and the
// integral increment is below 1/32, so every increment stays below 9 against a per-cycle bound of 15625.
TEMPLATE_TEST_CASE("velocity_form with an idle rate limiter emits the unlimited increments bit for bit",
    "[pid][siso][velocity-form][rate-limit]", float, double)
{
    using S = TestType;
    using LimitedPid = ctrlpp::pid<S, 1, ctrlpp::velocity_form, ctrlpp::rate_limit>;
    using FreePid = ctrlpp::pid<S, 1, ctrlpp::velocity_form>;
    constexpr std::array<double, 8> setpoints{1.0, 1.0, 1.5, 1.5, 2.0, 0.5, 0.5, 1.0};
    constexpr std::array<double, 8> measurements{0.0, 0.125, 0.375, 0.625, 1.0, 1.125, 0.875, 0.75};
    const auto gains = [](auto& cfg) { cfg.kp = vs(S{1}); cfg.ki = vs(S{1} / S{2}); cfg.kd = vs(S{1} / S{8}); };
    typename LimitedPid::config_type lcfg{};
    typename FreePid::config_type fcfg{};
    gains(lcfg);
    gains(fcfg);
    lcfg.template policy<ctrlpp::rate_limit>().rate_max = {static_cast<S>(1e6)};
    auto limited = ctrlpp::test::constructed(LimitedPid::create(lcfg));
    auto unlimited = ctrlpp::test::constructed(FreePid::create(fcfg));
    const S dt = S{1} / S{64};
    for(std::size_t k = 0; k < setpoints.size(); ++k)
    {
        CAPTURE(k);
        const auto sp = vs(static_cast<S>(setpoints[k]));
        const auto meas = vs(static_cast<S>(measurements[k]));
        CHECK(ctrlpp::test::commanded(limited.compute(sp, meas, dt))[0] == ctrlpp::test::commanded(unlimited.compute(sp, meas, dt))[0]);
    }
}
