#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <type_traits>
#include <initializer_list>

using ctrlpp::online_planner_substitution_reason;

namespace
{

template <typename Scalar>
struct limit_set
{
    Scalar v_max;
    Scalar a_max;
    Scalar j_max;
};

// The last three are limit sets on which a cruise integrated forward through its ramp rounds above v_max.
template <typename Scalar>
auto retarget_limit_sets() -> std::array<limit_set<Scalar>, 8>
{
    return {{
        {Scalar{1}, Scalar{1}, Scalar{1}},
        {Scalar{5}, Scalar{10}, Scalar{50}},
        {Scalar{2}, Scalar{3}, Scalar{7}},
        {static_cast<Scalar>(0.7), static_cast<Scalar>(1.3), static_cast<Scalar>(2.9)},
        {Scalar{10}, Scalar{100}, Scalar{1000}},
        {Scalar{1}, Scalar{5}, Scalar{50}},
        {Scalar{40}, Scalar{20}, Scalar{100}},
        {static_cast<Scalar>(1e-3), static_cast<Scalar>(1e-2), static_cast<Scalar>(1e-1)},
    }};
}

template <typename Scalar>
auto make_planner(limit_set<Scalar> const& limits) -> ctrlpp::online_planner_3rd<Scalar>
{
    auto created = ctrlpp::online_planner_3rd<Scalar>::create(
        {.v_max = limits.v_max, .a_max = limits.a_max, .j_max = limits.j_max});
    REQUIRE(created.has_value());
    return *std::move(created);
}

template <typename Scalar>
struct ramp_times
{
    Scalar T_j;
    Scalar T_a;
};

// @cite biagiotti2009 -- Sec. 3.4.3, rest-to-rest ramp with and without a_max reached
template <typename Scalar>
auto ramp_of(limit_set<Scalar> const& limits) -> ramp_times<Scalar>
{
    if (limits.v_max * limits.j_max >= limits.a_max * limits.a_max) {
        auto const T_j = limits.a_max / limits.j_max;
        return {T_j, T_j + limits.v_max / limits.a_max};
    }
    auto const T_j = std::sqrt(limits.v_max / limits.j_max);
    return {T_j, Scalar{2} * T_j};
}

enum class retarget_window : std::uint8_t
{
    cruise,
    acceleration_ramp,
    acceleration_nulling,
};

constexpr std::array<retarget_window, 3> all_windows{
    retarget_window::cruise, retarget_window::acceleration_ramp, retarget_window::acceleration_nulling};

template <typename Scalar>
auto window_bounds(retarget_window window, ramp_times<Scalar> ramp, Scalar dt) -> std::array<Scalar, 2>
{
    switch (window) {
    case retarget_window::cruise:
        return {Scalar{1.5} * ramp.T_a, Scalar{1.5} * ramp.T_a + Scalar{40} * dt};
    case retarget_window::acceleration_ramp:
        return {dt, ramp.T_a - dt};
    case retarget_window::acceleration_nulling:
        return {ramp.T_a - ramp.T_j, ramp.T_a - dt};
    }
    return {};
}

constexpr std::int32_t instants_per_window = 8;

template <typename Scalar>
auto window_instants(std::array<Scalar, 2> const& bounds, Scalar dt) -> std::array<Scalar, instants_per_window>
{
    auto const first = static_cast<std::int32_t>(std::ceil(bounds[0] / dt));
    auto const last = static_cast<std::int32_t>(std::floor(bounds[1] / dt));
    std::array<Scalar, instants_per_window> instants{};
    for (std::int32_t i = 0; i < instants_per_window; ++i) {
        auto const k = first + (last - first) * i / (instants_per_window - 1);
        instants[static_cast<std::size_t>(i)] = static_cast<Scalar>(k) * dt;
    }
    return instants;
}

struct plan_scan
{
    std::int32_t velocity_violations{};
    std::int32_t acceleration_violations{};
    std::int32_t jerk_violations{};
};

template <typename Scalar>
void record_limits(ctrlpp::trajectory_point<Scalar, 1> const& point, limit_set<Scalar> const& limits, plan_scan& scan)
{
    scan.velocity_violations += (std::abs(point.velocity(0)) > limits.v_max) ? 1 : 0;
    scan.acceleration_violations += (std::abs(point.acceleration(0)) > limits.a_max) ? 1 : 0;
}

// Elapsed time reaches a sample through at most 13 roundings: the reference subtraction, ten full-phase sums, the
// in-phase subtraction and, in a pinned phase, the time-to-end subtraction. Every phase group ends pinned, so
// acceleration carries at most 4 since the last pin: two phases integrated forward at a product and a sum each.
template <typename Scalar>
auto jerk_allowance(Scalar time_scale, limit_set<Scalar> const& limits) -> Scalar
{
    constexpr std::int32_t time_rounding_ops = 13;
    constexpr std::int32_t acceleration_rounding_ops = 4;
    auto const eps = std::numeric_limits<Scalar>::epsilon();
    return static_cast<Scalar>(time_rounding_ops) * eps * limits.j_max * time_scale
           + static_cast<Scalar>(acceleration_rounding_ops) * eps * limits.a_max;
}

template <typename Scalar>
void scan_dense(ctrlpp::online_planner_3rd<Scalar> const& planner, Scalar t_begin, Scalar dt, std::int32_t steps,
                limit_set<Scalar> const& limits, plan_scan& scan)
{
    auto const duration = planner.diagnostics().planned_duration;
    Scalar t_prev{};
    Scalar a_prev{};
    for (std::int32_t k = 0; k <= steps; ++k) {
        auto const t = t_begin + static_cast<Scalar>(k) * dt;
        auto const point = planner.sample(t);
        record_limits(point, limits, scan);
        auto const a = point.acceleration(0);
        auto const allowance = Scalar{2} * jerk_allowance(std::max(std::abs(t), duration), limits);
        if (k > 0 && std::abs(a - a_prev) > limits.j_max * (t - t_prev) + allowance) {
            ++scan.jerk_violations;
        }
        t_prev = t;
        a_prev = a;
    }
}

// Dense over the first and last 20 s of the plan, where every phase but the cruise lives, and every second between.
template <typename Scalar>
auto scan_plan(ctrlpp::online_planner_3rd<Scalar> const& planner, Scalar t0, limit_set<Scalar> const& limits, Scalar dt)
    -> plan_scan
{
    constexpr std::int32_t dense_seconds = 20;
    auto const duration = planner.diagnostics().planned_duration;
    auto const steps = static_cast<std::int32_t>(std::lround(static_cast<Scalar>(dense_seconds) / dt));
    plan_scan scan{};
    scan_dense(planner, t0, dt, steps, limits, scan);
    for (std::int32_t k = dense_seconds; static_cast<Scalar>(k + dense_seconds) < duration; ++k) {
        record_limits(planner.sample(t0 + static_cast<Scalar>(k)), limits, scan);
    }
    auto const tail = std::max(duration - static_cast<Scalar>(dense_seconds), Scalar{0});
    scan_dense(planner, t0 + tail, dt, steps, limits, scan);
    return scan;
}

template <typename Scalar>
struct sampled_planner
{
    ctrlpp::online_planner_3rd<Scalar> planner;
    ctrlpp::trajectory_point<Scalar, 1> point;
};

// A planner at rest at 0, commanded toward 1000 v_max, sampled once at t.
template <typename Scalar>
auto cruise_planner_sampled_at(limit_set<Scalar> const& limits, Scalar t) -> sampled_planner<Scalar>
{
    auto planner = make_planner(limits);
    planner.reset(Scalar{0});
    REQUIRE(planner.update(Scalar{1000} * limits.v_max).has_value());
    auto const point = planner.sample(t);
    return {std::move(planner), point};
}

constexpr double grid_dt = 1e-3;

template <typename Scalar>
struct sample_peaks
{
    Scalar velocity{};
    Scalar acceleration{};
    Scalar position{-std::numeric_limits<Scalar>::infinity()};
};

template <typename Scalar>
void record_peaks(ctrlpp::trajectory_point<Scalar, 1> const& point, sample_peaks<Scalar>& peaks)
{
    peaks.velocity = std::max(peaks.velocity, std::abs(point.velocity(0)));
    peaks.acceleration = std::max(peaks.acceleration, std::abs(point.acceleration(0)));
    peaks.position = std::max(peaks.position, point.position(0));
}

// Samples every representable instant from `before` steps below `center` through `after` steps above it.
template <typename Scalar>
auto scan_representable(ctrlpp::online_planner_3rd<Scalar> const& planner, Scalar center, std::int32_t before,
                        std::int32_t after) -> sample_peaks<Scalar>
{
    auto const infinity = std::numeric_limits<Scalar>::infinity();
    auto t = center;
    for (std::int32_t i = 0; i < before; ++i) {
        t = std::nextafter(t, -infinity);
    }
    sample_peaks<Scalar> peaks{};
    for (std::int32_t i = 0; i <= before + after; ++i) {
        record_peaks(planner.sample(t), peaks);
        t = std::nextafter(t, infinity);
    }
    return peaks;
}

constexpr std::int32_t representable_window = 64;

// So large in each type that the phase nulling a_max on unit limits moves the velocity and position by less than
// their resolution.
template <typename Scalar>
auto unresolvable_jerk() -> Scalar
{
    return std::is_same_v<Scalar, float> ? static_cast<Scalar>(1e8) : static_cast<Scalar>(1e16);
}

// j_max times the rounded quotient a_max / j_max rounds above a_max on this set. The double j_max is the fuzz input's
// 1e6 scaled by 2^-20, which scales the quotient and the product exactly and brings T_j to the sample grid's scale.
template <typename Scalar>
auto round_up_pair() -> std::array<Scalar, 2>
{
    if constexpr (std::is_same_v<Scalar, float>) {
        return {static_cast<Scalar>(0.1), Scalar{3}};
    } else {
        return {0x1.fbfbfbdbfbfbfp-4, 0x1.e848p-1};
    }
}

// The round-up pair with v_max = 4 a_max^2 / j_max, so every ramp and brake reaches a_max and holds it.
template <typename Scalar>
auto round_up_limits() -> limit_set<Scalar>
{
    auto const [a_max, j_max] = round_up_pair<Scalar>();
    REQUIRE(j_max * (a_max / j_max) > a_max);
    return {Scalar{4} * a_max * a_max / j_max, a_max, j_max};
}

}

TEST_CASE("OnlinePlanner3rd retarget: a cruise that integrated above the limit now carries it exactly",
          "[trajectory][online_planner_3rd][retarget]")
{
    auto const sets = std::array<limit_set<double>, 3>{{{1.0, 5.0, 50.0}, {40.0, 20.0, 100.0}, {1e-3, 1e-2, 1e-1}}};
    for (auto const& limits : sets) {
        CAPTURE(limits.v_max, limits.a_max, limits.j_max);
        auto planner = make_planner(limits);
        planner.reset(0.0);
        REQUIRE(planner.update(1000.0 * limits.v_max).has_value());

        auto const start = 1.5 * ramp_of(limits).T_a;
        double t = start;
        for (std::int32_t k = 0; k <= 40; ++k) {
            t = start + static_cast<double>(k) * grid_dt;
            auto const point = planner.sample(t);
            CAPTURE(k, t);
            REQUIRE(point.velocity(0) == limits.v_max);
            REQUIRE(point.acceleration(0) == 0.0);
        }

        REQUIRE(planner.update(2000.0 * limits.v_max).has_value());
        REQUIRE(planner.diagnostics().substitution_reason == online_planner_substitution_reason::none);
        auto const scan = scan_plan(planner, t, limits, grid_dt);
        CHECK(scan.velocity_violations == 0);
        CHECK(scan.acceleration_violations == 0);
        REQUIRE(planner.sample(t + 500.0).velocity(0) == limits.v_max);
    }
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: same-direction retargets from cruise, ramp and nulling stay within the "
                   "limits and never brake",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const dt = static_cast<Scalar>(grid_dt);
    std::int32_t retargets = 0;
    for (auto const& limits : retarget_limit_sets<Scalar>()) {
        for (auto const window : all_windows) {
            for (auto const t : window_instants(window_bounds(window, ramp_of(limits), dt), dt)) {
                CAPTURE(limits.v_max, limits.a_max, limits.j_max, static_cast<std::int32_t>(window), t);
                auto sampled = cruise_planner_sampled_at(limits, t);
                REQUIRE(sampled.planner.update(Scalar{2000} * limits.v_max).has_value());
                CHECK(sampled.planner.diagnostics().substitution_reason == online_planner_substitution_reason::none);
                auto const scan = scan_plan(sampled.planner, t, limits, dt);
                CHECK(scan.velocity_violations == 0);
                CHECK(scan.acceleration_violations == 0);
                ++retargets;
            }
        }
    }
    REQUIRE(retargets == 192);
}

TEST_CASE("OnlinePlanner3rd retarget: finite-difference jerk stays within j_max plus counted rounding",
          "[trajectory][online_planner_3rd][retarget]")
{
    std::int32_t retargets = 0;
    for (auto const& limits : retarget_limit_sets<double>()) {
        for (auto const window : all_windows) {
            for (auto const t : window_instants(window_bounds(window, ramp_of(limits), grid_dt), grid_dt)) {
                CAPTURE(limits.v_max, limits.a_max, limits.j_max, static_cast<std::int32_t>(window), t);
                auto sampled = cruise_planner_sampled_at(limits, t);
                REQUIRE(sampled.planner.update(2000.0 * limits.v_max).has_value());
                CHECK(scan_plan(sampled.planner, t, limits, grid_dt).jerk_violations == 0);
                ++retargets;
            }
        }
    }
    REQUIRE(retargets == 192);
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: a reversing retarget brakes and stays within the limits",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const dt = static_cast<Scalar>(grid_dt);
    std::int32_t retargets = 0;
    for (auto const& limits : retarget_limit_sets<Scalar>()) {
        for (auto const window : {retarget_window::cruise, retarget_window::acceleration_ramp}) {
            for (auto const t : window_instants(window_bounds(window, ramp_of(limits), dt), dt)) {
                CAPTURE(limits.v_max, limits.a_max, limits.j_max, static_cast<std::int32_t>(window), t);
                auto sampled = cruise_planner_sampled_at(limits, t);
                auto const target = sampled.point.position(0) - Scalar{1000} * limits.v_max;
                REQUIRE(sampled.planner.update(target).has_value());
                CHECK(sampled.planner.diagnostics().substitution_reason
                      == online_planner_substitution_reason::reversal_or_overshoot);
                auto const scan = scan_plan(sampled.planner, t, limits, dt);
                CHECK(scan.velocity_violations == 0);
                CHECK(scan.acceleration_violations == 0);
                ++retargets;
            }
        }
    }
    REQUIRE(retargets == 128);
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: a replan from a jerk-unresolvable acceleration keeps the limits and stops "
                   "at the target",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const limits = limit_set<Scalar>{Scalar{1}, Scalar{1}, unresolvable_jerk<Scalar>()};
    auto const target = Scalar{10};
    auto const t0 = static_cast<Scalar>(0.5);
    auto planner = make_planner(limits);
    planner.reset(Scalar{0});
    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.sample(t0).acceleration(0) > Scalar{0});
    REQUIRE(planner.update(target).has_value());

    auto const duration = planner.diagnostics().planned_duration;
    auto const dt = static_cast<Scalar>(grid_dt);
    auto const steps = static_cast<std::int32_t>(std::floor(duration / dt));
    auto peaks = scan_representable(planner, t0 + duration, representable_window, 0);
    for (std::int32_t k = 0; k <= steps; ++k) {
        record_peaks(planner.sample(t0 + static_cast<Scalar>(k) * dt), peaks);
    }
    CHECK(peaks.acceleration <= limits.a_max);
    CHECK(peaks.velocity <= limits.v_max);
    CHECK(peaks.position <= target);
}

TEST_CASE("OnlinePlanner3rd retarget: a replan from a speed below the speed floor keeps the velocity continuous",
          "[trajectory][online_planner_3rd][retarget]")
{
    // The decoded fuzz input that exposed a velocity step at the first pinned group end after the replan.
    constexpr double pos = 0x1.d3d2a2a2ad701p-924;
    constexpr double target = 0x1.ec8bf3d3d3dfcp-39;
    constexpr double dt = 0x1.0c6f7a0b5ed8dp-20;
    constexpr double t_start = 0x1.0624dd2f1a9fcp-9;
    constexpr double retarget_target = 0x1.f2a2a2a2a2a2ap-348;
    auto const limits = limit_set<double>{0x1.e848p+19, 0x1.0624dd2f1a9fcp-10, 0x1.0624dd2f1a9fcp-10};

    auto planner = make_planner(limits);
    planner.reset(pos);
    REQUIRE(planner.update(target).has_value());

    // A sampled speed carries at most 12 roundings since the last pinned group end, two forward phases of six, and
    // a sample instant 3: the product and sum forming it and the planner's reference subtraction.
    constexpr double velocity_rounding_ops = 12.0;
    constexpr double time_rounding_ops = 3.0;
    auto const eps = std::numeric_limits<double>::epsilon();
    double t_prev = t_start + dt;
    double v_prev = planner.sample(t_prev).velocity(0);
    REQUIRE(planner.update(retarget_target).has_value());
    double worst_excess = -limits.a_max;
    for (std::int32_t i = 2; i <= 2000; ++i) {
        auto const t = t_start + static_cast<double>(i) * dt;
        auto const v = planner.sample(t).velocity(0);
        auto const allowance = eps
                               * (velocity_rounding_ops * (std::abs(v) + std::abs(v_prev))
                                  + time_rounding_ops * limits.a_max * (std::abs(t) + std::abs(t_prev)))
                               / dt;
        worst_excess = std::max(worst_excess, std::abs(v - v_prev) / dt - limits.a_max - allowance);
        t_prev = t;
        v_prev = v;
    }
    CHECK(worst_excess <= 0.0);
}

TEST_CASE("OnlinePlanner3rd retarget: the decoded rest-to-rest fuzz input on a round-up set holds no acceleration "
          "above a_max",
          "[trajectory][online_planner_3rd][retarget]")
{
    // The decoded fuzz input whose ramp held a_max plus one unit in the last place through its constant phase.
    constexpr double pos = 0x0.00000ffffff7ep-1022;
    constexpr double target = 0x1.00046p-766;
    constexpr double dt = 0x1.0c6f7a0b5ed8dp-20;
    constexpr double t_start = -0x1.0624dd2f1a9fcp-9;
    constexpr double retarget_target = -0x1.0c1c1000096p+17;
    auto const limits = limit_set<double>{0x1.0624dd2f1a9fcp-10, 0x1.fbfbfbdbfbfbfp-4, 0x1.e848p+19};
    REQUIRE(limits.j_max * (limits.a_max / limits.j_max) > limits.a_max);

    auto planner = make_planner(limits);
    planner.reset(pos);
    REQUIRE(planner.update(target).has_value());
    auto const t0 = t_start + dt;
    REQUIRE(planner.sample(t0).velocity(0) == 0.0);
    REQUIRE(planner.update(retarget_target).has_value());

    // Only scan_plan's two dense windows: the cruise between them lasts about 1.4e8 s, where sampling it once a
    // second costs far more than the rest of the suite.
    auto const steps = static_cast<std::int32_t>(std::lround(20.0 / grid_dt));
    plan_scan scan{};
    scan_dense(planner, t0, grid_dt, steps, limits, scan);
    scan_dense(planner, t0 + planner.diagnostics().planned_duration - 20.0, grid_dt, steps, limits, scan);
    CHECK(scan.acceleration_violations == 0);
    CHECK(scan.velocity_violations == 0);
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: a rest-to-rest plan on a round-up set holds no acceleration above a_max",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const limits = round_up_limits<Scalar>();
    auto const dt = static_cast<Scalar>(grid_dt);
    auto planner = make_planner(limits);
    planner.reset(Scalar{0});
    REQUIRE(planner.update(-Scalar{1000} * limits.v_max).has_value());
    auto const scan = scan_plan(planner, Scalar{0}, limits, dt);
    CHECK(scan.acceleration_violations == 0);
    CHECK(scan.velocity_violations == 0);
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: a brake at the speed where its two shapes meet holds no acceleration "
                   "above a_max",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const limits = limit_set<Scalar>{static_cast<Scalar>(1e-3), static_cast<Scalar>(1e-2), static_cast<Scalar>(1e-1)};
    auto const dt = static_cast<Scalar>(grid_dt);
    for (std::int32_t const k : {107, 112, 117}) {
        auto const t = static_cast<Scalar>(k) * dt;
        CAPTURE(k, t);
        auto sampled = cruise_planner_sampled_at(limits, t);
        REQUIRE(sampled.planner.update(sampled.point.position(0) - Scalar{1000} * limits.v_max).has_value());
        CHECK(sampled.planner.diagnostics().substitution_reason
              == online_planner_substitution_reason::reversal_or_overshoot);
        auto const scan = scan_plan(sampled.planner, t, limits, dt);
        CHECK(scan.acceleration_violations == 0);
        CHECK(scan.velocity_violations == 0);
        auto const brake_peak = t + std::abs(sampled.point.acceleration(0)) / limits.j_max
                                + std::sqrt(limits.v_max / limits.j_max);
        auto const peaks = scan_representable(sampled.planner, brake_peak, representable_window, representable_window);
        CHECK(peaks.acceleration <= limits.a_max);
    }
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: a brake that reaches a_max on a round-up set holds no acceleration "
                   "above it",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const limits = round_up_limits<Scalar>();
    auto const dt = static_cast<Scalar>(grid_dt);
    auto const t = Scalar{1.5} * ramp_of(limits).T_a;
    auto sampled = cruise_planner_sampled_at(limits, t);
    REQUIRE(sampled.planner.update(sampled.point.position(0) - Scalar{1000} * limits.v_max).has_value());
    REQUIRE(sampled.planner.diagnostics().substitution_reason
            == online_planner_substitution_reason::reversal_or_overshoot);
    auto const scan = scan_plan(sampled.planner, t, limits, dt);
    CHECK(scan.acceleration_violations == 0);
    CHECK(scan.velocity_violations == 0);
}

TEMPLATE_TEST_CASE("OnlinePlanner3rd retarget: nulling an acceleration held at a_max on a round-up set stays within it",
                   "[trajectory][online_planner_3rd][retarget]", float, double)
{
    using Scalar = TestType;
    auto const limits = round_up_limits<Scalar>();
    auto const ramp = ramp_of(limits);
    auto const t = ramp.T_j + (ramp.T_a - Scalar{2} * ramp.T_j) / Scalar{2};
    auto sampled = cruise_planner_sampled_at(limits, t);
    CHECK(std::abs(sampled.point.acceleration(0)) <= limits.a_max);
    REQUIRE(sampled.planner.update(Scalar{1000} * limits.v_max).has_value());
    auto const peaks = scan_representable(sampled.planner, t, 0, representable_window);
    CHECK(peaks.acceleration <= limits.a_max);
}
