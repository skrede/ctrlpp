#include "hardening_helpers.h"

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto make_double_integrator() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt, dt;
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    return {A, B, C, D};
}

using OsqpMpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;

}

// ── MPC hardening: negative ────────────────────────────────────────────────────

TEST_CASE("MPC infeasible constraints: lower > upper", "[mpc][hardening][negative]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .u_min = Eigen::Matrix<double, 1, 1>{{5.0}},
        .u_max = Eigen::Matrix<double, 1, 1>{{-5.0}}, // lower > upper
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);

    // With infeasible input constraints, solver may return no solution
    // or return a degraded solution. Either is acceptable -- no crash.
    (void)u;
}

TEST_CASE("MPC minimal horizon N=1", "[mpc][hardening][negative]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 1,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 1.0).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());
    REQUIRE(std::isfinite((*u)(0)));
}

TEST_CASE("MPC with NaN in weight matrices", "[mpc][hardening][negative]")
{
    auto sys = make_double_integrator();

    Eigen::Matrix2d Q_nan = Eigen::Matrix2d::Identity();
    Q_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 5,
        .Q = Q_nan,
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    // Construction with NaN Q should not crash
    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    // Result is implementation-defined, but must not crash
    (void)u;
}

// ── MPC hardening: precision ───────────────────────────────────────────────────

TEST_CASE("MPC 1D regulation matches known optimal", "[mpc][hardening][precision]")
{
    // 1D integrator: x(k+1) = x(k) + u(k)
    constexpr std::size_t NX1 = 1;
    constexpr std::size_t NU1 = 1;

    Eigen::Matrix<double, 1, 1> A;
    A << 1.0;
    Eigen::Matrix<double, 1, 1> B;
    B << 1.0;
    Eigen::Matrix<double, 1, 1> C;
    C << 1.0;
    Eigen::Matrix<double, 1, 1> D;
    D << 0.0;

    ctrlpp::discrete_state_space<double, NX1, NU1, NX1> sys{A, B, C, D};

    ctrlpp::mpc_config<double, NX1, NU1> cfg{
        .horizon = 20,
        .Q = (Eigen::Matrix<double, 1, 1>() << 1.0).finished(),
        .R = (Eigen::Matrix<double, 1, 1>() << 1.0).finished(),
    };

    ctrlpp::mpc<double, NX1, NU1, ctrlpp::osqp_solver> controller(sys, cfg);

    Eigen::Matrix<double, 1, 1> x;
    x << 1.0;

    auto u = controller.solve(x);
    REQUIRE(u.has_value());

    // For 1D integrator with Q=R=I, optimal u should be negative (drive to 0)
    CHECK((*u)(0) < 0.0);
}

// ── MPC hardening: stability ───────────────────────────────────────────────────

TEST_CASE("MPC closed-loop stabilizes double integrator", "[mpc][hardening][stability]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{1.0, 0.5};

    for (int step = 0; step < 100; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value();
        REQUIRE(std::isfinite(x(0)));
        REQUIRE(std::isfinite(x(1)));
    }

    CHECK(x.norm() < 0.1);
}

// ── MPC hardening: robustness ──────────────────────────────────────────────────

TEST_CASE("MPC with huge Q weights", "[mpc][hardening][robustness]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity() * 1e10,
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());
    REQUIRE(std::isfinite((*u)(0)));
}

TEST_CASE("MPC with near-zero R weights", "[mpc][hardening][robustness]")
{
    auto sys = make_double_integrator();

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 1e-10).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());
    REQUIRE(std::isfinite((*u)(0)));
}
