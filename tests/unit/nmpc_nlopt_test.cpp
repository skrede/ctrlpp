#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto double_integrator = [](const Eigen::Vector2d& x,
                            const Eigen::Matrix<double, 1, 1>& u)
    -> Eigen::Vector2d {
    return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)};
};

auto pendulum = [](const Eigen::Vector2d& x,
                   const Eigen::Matrix<double, 1, 1>& u)
    -> Eigen::Vector2d {
    constexpr double pdt = 0.05;
    constexpr double g = 9.81;
    constexpr double l = 1.0;
    double theta = x(0);
    double omega = x(1);
    double alpha = -g / l * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + pdt * omega, omega + pdt * alpha};
};

auto make_config(int horizon = 10) -> ctrlpp::nmpc_config<double, NX, NU> {
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };
}

using NloptSolver = ctrlpp::nlopt_solver<double>;
using NmpcDI = ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator)>;
using NmpcPend = ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(pendulum)>;

}

TEST_CASE("nlopt_solver satisfies nlp_solver concept") {
    static_assert(ctrlpp::nlp_solver<NloptSolver>,
        "nlopt_solver<double> must satisfy nlp_solver concept");
}

TEST_CASE("nmpc nlopt regulation", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    double initial_norm = x.norm();

    for (int step = 0; step < 50; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    REQUIRE(x.norm() < 0.1 * initial_norm);
}

TEST_CASE("nmpc nlopt setpoint tracking", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.0, 0.0};
    Eigen::Vector2d x_ref{2.0, 0.0};

    for (int step = 0; step < 80; ++step) {
        auto u = controller.solve(x, x_ref);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    REQUIRE((x - x_ref).norm() < 0.5);
}

TEST_CASE("nmpc nlopt input box constraints", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();
    config.u_min = Eigen::Matrix<double, 1, 1>{-0.5};
    config.u_max = Eigen::Matrix<double, 1, 1>{0.5};

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{5.0, 0.0};

    for (int step = 0; step < 30; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        CHECK((*u)(0) >= -0.5 - 1e-6);
        CHECK((*u)(0) <= 0.5 + 1e-6);
        x = double_integrator(x, *u);
    }
}

TEST_CASE("nmpc nlopt state box constraints", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();
    config.x_min = Eigen::Vector2d{-2.0, -2.0};
    config.x_max = Eigen::Vector2d{2.0, 2.0};

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.5, 0.0};

    for (int step = 0; step < 30; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);

        auto [states, inputs] = controller.trajectory();
        for (const auto& s : states) {
            CHECK(s(0) >= -2.0 - 1e-4);
            CHECK(s(0) <= 2.0 + 1e-4);
            CHECK(s(1) >= -2.0 - 1e-4);
            CHECK(s(1) <= 2.0 + 1e-4);
        }
    }
}

TEST_CASE("nmpc nlopt rate constraints", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();
    config.du_max = Eigen::Matrix<double, 1, 1>{0.1};

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{2.0, 0.0};
    double u_prev = 0.0;

    for (int step = 0; step < 20; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        double du = std::abs((*u)(0) - u_prev);
        CHECK(du <= 0.1 + 1e-4);

        u_prev = (*u)(0);
        x = double_integrator(x, *u);
    }
}

TEST_CASE("nmpc nlopt warm-start benefit", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    // First solve (cold start)
    auto u1 = controller.solve(x);
    REQUIRE(u1.has_value());
    auto diag1 = controller.diagnostics();

    // Step forward
    x = double_integrator(x, *u1);

    // Second solve (warm start from shifted solution)
    auto u2 = controller.solve(x);
    REQUIRE(u2.has_value());
    auto diag2 = controller.diagnostics();

    // Warm start should use fewer or equal evaluations
    CHECK(diag2.iterations <= diag1.iterations);
}

TEST_CASE("nmpc nlopt custom cost", "[nmpc][nlopt]") {
    // Default quadratic cost
    auto default_config = make_config(10);
    default_config.Q = Eigen::Matrix2d::Identity();
    default_config.R = Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI default_ctrl{double_integrator, default_config};

    Eigen::Vector2d x0{1.0, 1.0};
    auto u_default = default_ctrl.solve(x0);
    REQUIRE(u_default.has_value());

    // Custom cost: heavily penalize position, ignore velocity
    auto custom_config = make_config(10);
    custom_config.stage_cost = [](const Eigen::Vector2d& x,
                                  const Eigen::Matrix<double, 1, 1>& u) -> double {
        return 100.0 * x(0) * x(0) + 0.01 * u(0) * u(0);
    };
    custom_config.terminal_cost = [](const Eigen::Vector2d& x) -> double {
        return 100.0 * x(0) * x(0);
    };

    NmpcDI custom_ctrl{double_integrator, custom_config};
    auto u_custom = custom_ctrl.solve(x0);
    REQUIRE(u_custom.has_value());

    // Custom cost should produce different control due to different weighting
    CHECK(std::abs((*u_default)(0) - (*u_custom)(0)) > 1e-3);
}

TEST_CASE("nmpc nlopt trajectory tracking", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.0, 0.0};
    double max_error = 0.0;

    for (int step = 0; step < 50; ++step) {
        double t = step * dt;

        // Sinusoidal reference trajectory
        std::vector<Eigen::Vector2d> refs;
        refs.reserve(11);
        for (int k = 0; k <= 10; ++k) {
            double tk = t + k * dt;
            refs.push_back(Eigen::Vector2d{
                std::sin(0.5 * tk),
                0.5 * std::cos(0.5 * tk)
            });
        }

        auto u = controller.solve(x,
            std::span<const Eigen::Vector2d>{refs});
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);

        double error = (x - refs[1]).norm();
        max_error = std::max(max_error, error);
    }

    // Tracking error should be bounded
    CHECK(max_error < 2.0);
}

TEST_CASE("nmpc nlopt pendulum regulation", "[nmpc][nlopt]") {
    auto config = make_config(15);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.01 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcPend controller{pendulum, config};

    Eigen::Vector2d x{0.5, 0.0};
    double initial_norm = x.norm();

    for (int step = 0; step < 80; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = pendulum(x, *u);
    }

    REQUIRE(x.norm() < 0.2 * initial_norm);
}

TEST_CASE("nmpc nlopt constraint satisfaction closed-loop", "[nmpc][nlopt]") {
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();
    config.u_min = Eigen::Matrix<double, 1, 1>{-1.0};
    config.u_max = Eigen::Matrix<double, 1, 1>{1.0};
    config.x_min = Eigen::Vector2d{-3.0, -3.0};
    config.x_max = Eigen::Vector2d{3.0, 3.0};

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{2.5, 0.5};
    constexpr double tol = 1e-4;

    for (int step = 0; step < 100; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        // Input constraint check
        CHECK((*u)(0) >= -1.0 - tol);
        CHECK((*u)(0) <= 1.0 + tol);

        x = double_integrator(x, *u);

        // State constraint check
        CHECK(x(0) >= -3.0 - tol);
        CHECK(x(0) <= 3.0 + tol);
        CHECK(x(1) >= -3.0 - tol);
        CHECK(x(1) <= 3.0 + tol);
    }
}
