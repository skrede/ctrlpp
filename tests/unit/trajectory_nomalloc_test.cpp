// Verify the steady-state trajectory hot paths do zero heap allocation on
// fixed-size templated inputs, using the belt-and-suspenders harness from
// nomalloc_harness.h: a throwing eigen_assert that survives -DNDEBUG plus a
// global allocation counter that catches heap traffic outside Eigen's own
// bookkeeping. The harness header must stay the first include of this file.
//
// Coverage: evaluate() for cubic/quintic/septic polynomial segments, the
// trapezoidal/double_s/modified_sin/modified_trap profiles, and the
// cubic_spline/smoothing_spline/bspline splines, plus online_planner_2nd/3rd
// update()+sample(), the time rescaling of both velocity profiles, and both
// multi-axis synchronization overloads. Construction (splines solve linear
// systems at build time by design) happens outside the armed window; only the
// steady-state evaluate, sample, retime or synchronize call is guarded.

#include "nomalloc_harness.h"

#include "ctrlpp/trajectory/synchronize.h"
#include "ctrlpp/trajectory/cubic_spline.h"
#include "ctrlpp/trajectory/smoothing_spline.h"
#include "ctrlpp/trajectory/cubic_trajectory.h"
#include "ctrlpp/trajectory/septic_trajectory.h"
#include "ctrlpp/trajectory/quintic_trajectory.h"
#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/bspline_trajectory.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/modified_sin_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"
#include "ctrlpp/trajectory/modified_trap_trajectory.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <span>
#include <array>
#include <vector>
#include <cstddef>
#include <utility>


namespace
{

using ctrlpp::Vector;

template <typename Fn>
std::size_t guarded_allocations(Fn&& fn)
{
    ctrlpp_test::scoped_no_malloc guard;
    std::forward<Fn>(fn)();
    return guard.allocations();
}

// create() is the only construction path on the two velocity profiles. It hands
// back a ctrlpp::expected holding the profile by value, so it allocates nothing
// either, but construction still happens outside the armed window here for the
// same reason every other construction in this file does: only the steady-state
// call is under test.
auto trapezoidal_axis_profile(ctrlpp::trapezoidal_trajectory<double>::config const& cfg)
    -> ctrlpp::trapezoidal_trajectory<double>
{
    auto created = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

auto double_s_axis_profile(ctrlpp::double_s_trajectory<double>::config const& cfg)
    -> ctrlpp::double_s_trajectory<double>
{
    auto created = ctrlpp::double_s_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

// Dense-sample a segment's evaluate() across [0, duration] inside the armed
// window and return the observed allocation count.
template <typename Segment>
std::size_t evaluate_sweep_allocations(const Segment& seg, double duration, int samples)
{
    return guarded_allocations([&] {
        for(int i = 0; i <= samples; ++i)
            seg.evaluate(duration * static_cast<double>(i) / static_cast<double>(samples));
    });
}

}


TEST_CASE("cubic_trajectory evaluate performs zero heap allocation",
          "[trajectory][cubic][hardening][nomalloc]")
{
    const Vector<double, 1> q0 = Vector<double, 1>::Zero();
    const Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    const Vector<double, 1> v0 = Vector<double, 1>::Zero();
    const Vector<double, 1> v1 = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 2.0).value();

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("quintic_trajectory evaluate performs zero heap allocation",
          "[trajectory][quintic][hardening][nomalloc]")
{
    const Vector<double, 1> q0 = Vector<double, 1>::Zero();
    const Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    const Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_quintic_trajectory(q0, q1, zero, zero, zero, zero, 2.0).value();

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("septic_trajectory evaluate performs zero heap allocation",
          "[trajectory][septic][hardening][nomalloc]")
{
    const Vector<double, 1> q0 = Vector<double, 1>::Zero();
    const Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    const Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_septic_trajectory(q0, q1, zero, zero, zero, zero, zero, zero, 2.0).value();

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("trapezoidal_trajectory evaluate performs zero heap allocation",
          "[trajectory][trapezoidal][hardening][nomalloc]")
{
    auto const seg = trapezoidal_axis_profile({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("double_s_trajectory evaluate performs zero heap allocation",
          "[trajectory][double_s][hardening][nomalloc]")
{
    auto const seg =
        double_s_axis_profile({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("modified_sin_trajectory evaluate performs zero heap allocation",
          "[trajectory][modified_sin][hardening][nomalloc]")
{
    ctrlpp::modified_sin_trajectory<double> seg({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("modified_trap_trajectory evaluate performs zero heap allocation",
          "[trajectory][modified_trap][hardening][nomalloc]")
{
    ctrlpp::modified_trap_trajectory<double> seg({.q0 = 0.0, .q1 = 10.0, .T = 2.0});

    seg.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(seg, seg.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("cubic_spline evaluate performs zero heap allocation",
          "[trajectory][cubic_spline][hardening][nomalloc]")
{
    auto created = ctrlpp::cubic_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.5, 2.0, 1.5},
    });
    REQUIRE(created.has_value());
    auto& spline = *created;

    spline.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(spline, 4.0, 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("smoothing_spline evaluate performs zero heap allocation",
          "[trajectory][smoothing_spline][hardening][nomalloc]")
{
    auto created = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.5, 2.0, 1.5},
        .mu = 0.5,
    });
    REQUIRE(created.has_value());
    auto& spline = *created;

    spline.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(spline, 4.0, 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("bspline_trajectory evaluate performs zero heap allocation",
          "[trajectory][bspline][hardening][nomalloc]")
{
    auto created = ctrlpp::bspline_trajectory<double, 3>::try_create({
        .control_points = {0.0, 1.0, 3.0, 2.0, 4.0},
    });
    REQUIRE(created.has_value());
    auto& spline = *created;

    spline.evaluate(0.0);

    std::size_t allocations = 0;
    allocations = evaluate_sweep_allocations(spline, spline.duration(), 256);
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
}

TEST_CASE("online_planner_2nd update and sample perform zero heap allocation",
          "[trajectory][online_planner_2nd][hardening][nomalloc]")
{
    auto created = ctrlpp::online_planner_2nd<double>::try_create({.v_max = 5.0, .a_max = 10.0});
    REQUIRE(created.has_value());
    auto& planner = *created;

    planner.update(10.0);
    planner.sample(0.5);

    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
        {
            planner.update((i % 2 == 0) ? 5.0 : -5.0);
            planner.sample(0.01 * static_cast<double>(i));
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
}

namespace
{

// Axis configurations reused by the retiming and synchronization cases below.
// The trapezoidal pair drives the closed-form retiming path, the double-S pair
// drives the rest-to-rest closed form and the bracketed solve respectively.
constexpr ctrlpp::trapezoidal_trajectory<double>::config trapezoidal_axis{
    .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .v0 = 0.5, .v1 = 0.25};

constexpr ctrlpp::double_s_trajectory<double>::config double_s_rest_to_rest{
    .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0, .v0 = 0.0, .v1 = 0.0};

constexpr ctrlpp::double_s_trajectory<double>::config double_s_boundary_velocities{
    .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0, .v0 = 1.0, .v1 = 0.5};

}

TEST_CASE("trapezoidal_trajectory rescale_to performs zero heap allocation",
          "[trajectory][trapezoidal][hardening][nomalloc]")
{
    // Warm-up on a throwaway copy, outside the armed window: the retiming
    // mutates the profile, so the guarded call has to run on a fresh one.
    auto warmup = trapezoidal_axis_profile(trapezoidal_axis);
    REQUIRE(warmup.rescale_to(warmup.duration() * 2.0).has_value());

    auto seg = trapezoidal_axis_profile(trapezoidal_axis);
    const double target = seg.duration() * 2.0;

    bool retimed = false;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] { retimed = seg.rescale_to(target).has_value(); });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
    REQUIRE(retimed);
}

TEST_CASE("double_s_trajectory rescale_to performs zero heap allocation on both paths",
          "[trajectory][double_s][hardening][nomalloc]")
{
    auto rest_warmup = double_s_axis_profile(double_s_rest_to_rest);
    REQUIRE(rest_warmup.rescale_to(rest_warmup.duration() * 2.0).has_value());
    auto solved_warmup = double_s_axis_profile(double_s_boundary_velocities);
    REQUIRE(solved_warmup.rescale_to(solved_warmup.duration() * 1.2).has_value());

    // Rest to rest: the scale follows in closed form from the two durations.
    auto rest = double_s_axis_profile(double_s_rest_to_rest);
    const double rest_target = rest.duration() * 2.0;

    bool rest_retimed = false;
    std::size_t rest_allocations = 0;
    rest_allocations = guarded_allocations([&] { rest_retimed = rest.rescale_to(rest_target).has_value(); });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(rest_allocations == 0);
    REQUIRE(rest_retimed);

    // Nonzero boundary velocities: the scale is found by halving a bracket, and
    // every candidate profile the search builds lives on the stack.
    auto solved = double_s_axis_profile(double_s_boundary_velocities);
    const double solved_target = solved.duration() * 1.2;

    bool solved_retimed = false;
    std::size_t solved_allocations = 0;
    solved_allocations =
        guarded_allocations([&] { solved_retimed = solved.rescale_to(solved_target).has_value(); });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(solved_allocations == 0);
    REQUIRE(solved_retimed);
}

TEST_CASE("synchronize performs zero heap allocation, variadic overload",
          "[trajectory][hardening][nomalloc]")
{
    auto warm1 = trapezoidal_axis_profile({.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0});
    auto warm2 = trapezoidal_axis_profile({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0});
    auto warm3 = double_s_axis_profile(double_s_rest_to_rest);
    REQUIRE(ctrlpp::synchronize(warm1, warm2).has_value());
    REQUIRE(ctrlpp::synchronize(warm3, warm3).has_value());

    auto ax1 = trapezoidal_axis_profile({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0});
    auto ax2 = trapezoidal_axis_profile({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0});
    auto ax3 = trapezoidal_axis_profile({.q0 = 0.0, .q1 = 1.0, .v_max = 5.0, .a_max = 10.0});

    bool synchronized = false;
    std::size_t allocations = 0;
    allocations =
        guarded_allocations([&] { synchronized = ctrlpp::synchronize(ax1, ax2, ax3).has_value(); });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
    REQUIRE(synchronized);
}

TEST_CASE("synchronize performs zero heap allocation over a fixed-size backing store",
          "[trajectory][hardening][nomalloc]")
{
    // The contiguous-view overload exists for exactly this shape: a statically
    // allocated run of axes, no owning container anywhere on the path.
    using axis = ctrlpp::trapezoidal_trajectory<double>;

    std::array<axis, 3> warmup{
        trapezoidal_axis_profile({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0}),
        trapezoidal_axis_profile({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0}),
        trapezoidal_axis_profile({.q0 = 0.0, .q1 = 1.0, .v_max = 5.0, .a_max = 10.0})};
    REQUIRE(ctrlpp::synchronize(std::span<axis>{warmup}).has_value());

    std::array<axis, 3> axes{
        trapezoidal_axis_profile({.q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 10.0}),
        trapezoidal_axis_profile({.q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0}),
        trapezoidal_axis_profile({.q0 = 0.0, .q1 = 1.0, .v_max = 5.0, .a_max = 10.0})};

    bool synchronized = false;
    std::size_t allocations = 0;
    allocations = guarded_allocations(
        [&] { synchronized = ctrlpp::synchronize(std::span<axis>{axes}).has_value(); });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(allocations == 0);
    REQUIRE(synchronized);
}

TEST_CASE("online_planner_3rd update and sample perform zero heap allocation",
          "[trajectory][online_planner_3rd][hardening][nomalloc]")
{
    auto created =
        ctrlpp::online_planner_3rd<double>::try_create({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});
    REQUIRE(created.has_value());
    auto& planner = *created;

    planner.update(10.0);
    planner.sample(0.5);

    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
        {
            planner.update((i % 2 == 0) ? 5.0 : -5.0);
            planner.sample(0.01 * static_cast<double>(i));
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
}
