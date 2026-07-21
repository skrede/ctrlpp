// Verify the steady-state control-tier hot paths do zero heap allocation on
// fixed-size templated inputs, using the belt-and-suspenders harness from
// nomalloc_harness.h: a throwing eigen_assert that survives -DNDEBUG plus a
// global allocation counter that catches heap traffic outside Eigen's own
// bookkeeping. The harness header must stay the first include of this file.
//
// Coverage: pid compute() in both position and velocity form (including a
// composed anti-windup + derivative-filter variant), plus the infinite-horizon
// LQR steady-state control laws lqr::compute and lqr_time_varying::compute, the
// fixed-size per-tick matrix-vector multiplies the STACK-RT-SAFETY-CONTRACT
// groups with dare/care. Construction and gain computation happen outside the
// armed window; only the steady-state call is guarded.

#include "nomalloc_harness.h"

#include "ctrlpp/pid.h"

#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <vector>
#include <cstddef>
#include <utility>


namespace
{

// Runs the callable inside an armed no-malloc window and returns the number of
// heap allocations it performed. The count is sampled before the guard is
// released and before any test macro runs, so framework-internal allocations
// cannot pollute it.
template <typename Fn>
std::size_t guarded_allocations(Fn&& fn)
{
    ctrlpp_test::scoped_no_malloc guard;
    std::forward<Fn>(fn)();
    return guard.allocations();
}

}


TEST_CASE("pid position-form compute performs zero heap allocation",
          "[pid][hardening][nomalloc]")
{
    using pid_type = ctrlpp::pid<double, 1>;

    pid_type::config_type cfg{};
    cfg.kp = ctrlpp::Vector<double, 1>::Constant(2.5);
    cfg.ki = ctrlpp::Vector<double, 1>::Constant(0.5);
    cfg.kd = ctrlpp::Vector<double, 1>::Constant(0.1);
    pid_type controller(cfg);

    const auto sp = ctrlpp::Vector<double, 1>::Constant(1.0);
    const auto meas = ctrlpp::Vector<double, 1>::Zero();
    constexpr double dt = 0.01;

    for(int i = 0; i < 8; ++i)
        controller.compute(sp, meas, dt);

    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            controller.compute(sp, meas, dt);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
}

TEST_CASE("pid velocity-form compute performs zero heap allocation",
          "[pid][hardening][nomalloc]")
{
    using pid_type = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;

    pid_type::config_type cfg{};
    cfg.kp = ctrlpp::Vector<double, 1>::Constant(2.0);
    cfg.ki = ctrlpp::Vector<double, 1>::Constant(1.0);
    cfg.kd = ctrlpp::Vector<double, 1>::Constant(0.05);
    pid_type controller(cfg);

    const auto sp = ctrlpp::Vector<double, 1>::Constant(1.0);
    const auto meas = ctrlpp::Vector<double, 1>::Zero();
    constexpr double dt = 0.01;

    for(int i = 0; i < 8; ++i)
        controller.compute(sp, meas, dt);

    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            controller.compute(sp, meas, dt);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
}

TEST_CASE("pid composed anti-windup and derivative-filter compute performs zero heap allocation",
          "[pid][hardening][nomalloc]")
{
    using pid_type =
        ctrlpp::pid<double, 3, ctrlpp::anti_windup<ctrlpp::back_calc>, ctrlpp::deriv_filter>;

    pid_type::config_type cfg{};
    cfg.kp = ctrlpp::Vector<double, 3>::Constant(2.0);
    cfg.ki = ctrlpp::Vector<double, 3>::Constant(0.5);
    cfg.kd = ctrlpp::Vector<double, 3>::Constant(0.1);
    cfg.output_min = ctrlpp::Vector<double, 3>::Constant(-1.0);
    cfg.output_max = ctrlpp::Vector<double, 3>::Constant(1.0);
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {1.0, 1.0, 1.0};
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0, 10.0, 10.0};
    pid_type controller(cfg);

    const auto sp = ctrlpp::Vector<double, 3>::Constant(5.0);
    const auto meas = ctrlpp::Vector<double, 3>::Zero();
    constexpr double dt = 0.01;

    for(int i = 0; i < 8; ++i)
        controller.compute(sp, meas, dt);

    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            controller.compute(sp, meas, dt);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
}

TEST_CASE("lqr steady-state control law performs zero heap allocation",
          "[lqr][hardening][nomalloc]")
{
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;

    Eigen::Matrix<double, NX, NX> A = Eigen::Matrix<double, NX, NX>::Identity();
    for(int i = 0; i + 1 < int(NX); ++i)
        A(i, i + 1) = 0.05;

    Eigen::Matrix<double, NX, NU> B = Eigen::Matrix<double, NX, NU>::Zero();
    B(1, 0) = 0.05;
    B(3, 1) = 0.05;

    Eigen::Matrix<double, NX, NX> Q = Eigen::Matrix<double, NX, NX>::Identity();
    Eigen::Matrix<double, NU, NU> R = 0.1 * Eigen::Matrix<double, NU, NU>::Identity();

    auto gain = ctrlpp::lqr_gain<double, NX, NU>(A, B, Q, R);
    REQUIRE(gain.has_value());

    ctrlpp::lqr<double, NX, NU> controller(*gain);

    ctrlpp::lqr<double, NX, NU>::state_type x;
    x << 1.0, -0.5, 0.25, -0.125;

    for(int i = 0; i < 8; ++i)
        controller.compute(x);

    std::size_t allocations = 0;
    ctrlpp::lqr<double, NX, NU>::input_type u = ctrlpp::lqr<double, NX, NU>::input_type::Zero();
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            u = controller.compute(x);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(u.allFinite());
}

TEST_CASE("lqr_time_varying steady-state control law performs zero heap allocation",
          "[lqr][hardening][nomalloc]")
{
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;

    using gain_type = ctrlpp::lqr_time_varying<double, NX, NU>::gain_type;

    Eigen::Matrix<double, NX, NX> A = Eigen::Matrix<double, NX, NX>::Identity();
    for(int i = 0; i + 1 < int(NX); ++i)
        A(i, i + 1) = 0.05;

    Eigen::Matrix<double, NX, NU> B = Eigen::Matrix<double, NX, NU>::Zero();
    B(1, 0) = 0.05;
    B(3, 1) = 0.05;

    Eigen::Matrix<double, NX, NX> Q = Eigen::Matrix<double, NX, NX>::Identity();
    Eigen::Matrix<double, NU, NU> R = 0.1 * Eigen::Matrix<double, NU, NU>::Identity();

    auto steady = ctrlpp::lqr_gain<double, NX, NU>(A, B, Q, R);
    REQUIRE(steady.has_value());

    constexpr std::size_t horizon = 16;
    std::vector<gain_type> gains(horizon, *steady);
    ctrlpp::lqr_time_varying<double, NX, NU> controller(std::move(gains));

    ctrlpp::lqr_time_varying<double, NX, NU>::state_type x;
    x << 1.0, -0.5, 0.25, -0.125;

    for(std::size_t k = 0; k < horizon; ++k)
        controller.compute(x, k);

    std::size_t allocations = 0;
    ctrlpp::lqr_time_varying<double, NX, NU>::input_type u =
        ctrlpp::lqr_time_varying<double, NX, NU>::input_type::Zero();
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            u = controller.compute(x, static_cast<std::size_t>(i) % horizon);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(u.allFinite());
}
