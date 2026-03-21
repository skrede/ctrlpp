#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/state_space.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>

namespace {

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto make_double_integrator() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    Eigen::Matrix2d A;
    A << 1.0, dt,
         0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt,
         dt;
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    return {A, B, C, D};
}

auto double_integrator_dynamics = [](const Eigen::Vector2d& x,
                                     const Eigen::Matrix<double, 1, 1>& u)
    -> Eigen::Vector2d {
    return Eigen::Vector2d{x(0) + dt * x(1) + 0.5 * dt * dt * u(0),
                           x(1) + dt * u(0)};
};

using OsqpMpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;
using NloptSolver = ctrlpp::nlopt_solver<double>;
using NmpcDI = ctrlpp::nmpc<double, NX, NU, NloptSolver,
                             decltype(double_integrator_dynamics)>;

}

TEST_CASE("linear mpc closed-loop convergence - double integrator",
          "[mpc][closedloop]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.01).finished(),
        .u_min = (Eigen::Matrix<double, 1, 1>() << -10.0).finished(),
        .u_max = (Eigen::Matrix<double, 1, 1>() << 10.0).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{5.0, 0.0};
    Eigen::Vector2d x_ref{0.0, 0.0};

    for (int step = 0; step < 200; ++step) {
        auto u = controller.solve(x, x_ref);
        REQUIRE(u.has_value());
        REQUIRE(std::isfinite((*u)(0)));
        x = sys.A * x + sys.B * u.value();
        REQUIRE(std::isfinite(x(0)));
        REQUIRE(std::isfinite(x(1)));
    }

    CHECK(x.norm() < 0.1);
}

TEST_CASE("nonlinear mpc closed-loop convergence - double integrator",
          "[nmpc][closedloop]")
{
    ctrlpp::nmpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
        .u_min = (Eigen::Matrix<double, 1, 1>() << -10.0).finished(),
        .u_max = (Eigen::Matrix<double, 1, 1>() << 10.0).finished(),
    };

    NmpcDI controller{double_integrator_dynamics, cfg};

    Eigen::Vector2d x{5.0, 0.0};

    for (int step = 0; step < 200; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator_dynamics(x, *u);
    }

    CHECK(x.norm() < 0.5);
}

TEST_CASE("linear mpc trajectory tracking with reference change",
          "[mpc][closedloop]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.01).finished(),
        .u_min = (Eigen::Matrix<double, 1, 1>() << -10.0).finished(),
        .u_max = (Eigen::Matrix<double, 1, 1>() << 10.0).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{0.0, 0.0};
    Eigen::Vector2d ref1{0.0, 0.0};
    Eigen::Vector2d ref2{3.0, 0.0};

    // Phase 1: regulate to origin for 100 steps
    for (int step = 0; step < 100; ++step) {
        auto u = controller.solve(x, ref1);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value();
    }

    CHECK(x.norm() < 0.1);

    // Phase 2: track new reference for 200 steps
    for (int step = 100; step < 300; ++step) {
        auto u = controller.solve(x, ref2);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value();
    }

    CHECK((x - ref2).norm() < 0.1);
}
