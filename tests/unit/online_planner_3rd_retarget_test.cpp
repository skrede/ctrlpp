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
