#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>

namespace {

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto make_double_integrator() -> ctrlpp::DiscreteStateSpace<double, NX, NU, NX> {
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

using OsqpMpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;

}

TEST_CASE("mpc with OSQP solver", "[mpc][osqp]") {
    auto sys = make_double_integrator();
    constexpr int N = 10;

    SECTION("regulation drives state toward origin") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x{1.0, 0.0};
        double prev_norm = x.norm();

        // Simulate closed-loop for 50 steps
        int monotonic_after = 5; // allow transient in first few steps
        for (int step = 0; step < 50; ++step) {
            auto u = controller.solve(x);
            REQUIRE(u.has_value());
            x = sys.A * x + sys.B * u.value();

            double norm = x.norm();
            if (step >= monotonic_after) {
                CHECK(norm <= prev_norm + 1e-6);
            }
            prev_norm = norm;
        }

        CHECK(x.norm() < 0.1);
    }

    SECTION("setpoint tracking converges to reference") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x{0.0, 0.0};
        Eigen::Vector2d x_ref{2.0, 0.0};

        // First solve should produce non-zero control pushing toward reference
        auto u0 = controller.solve(x, x_ref);
        REQUIRE(u0.has_value());
        CHECK(std::abs((*u0)(0)) > 1e-6);

        // Simulate closed-loop for 50 steps
        for (int step = 0; step < 50; ++step) {
            auto u = controller.solve(x, x_ref);
            REQUIRE(u.has_value());
            x = sys.A * x + sys.B * u.value();
        }

        CHECK((x - x_ref).norm() < 0.1);
    }

    SECTION("input bounds are respected") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << -0.5).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << 0.5).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x{5.0, 0.0}; // large initial state to push solver hard

        for (int step = 0; step < 20; ++step) {
            auto u = controller.solve(x);
            REQUIRE(u.has_value());
            CHECK((*u)(0) >= -0.5 - 1e-4);
            CHECK((*u)(0) <= 0.5 + 1e-4);
            x = sys.A * x + sys.B * u.value();
        }
    }

    SECTION("soft state constraints allow feasible solution from outside bounds") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .x_min = (Eigen::Vector2d() << -5.0, -5.0).finished(),
            .x_max = (Eigen::Vector2d() << 5.0, 5.0).finished(),
            .hard_state_constraints = false,
        };

        Eigen::Vector2d x0{10.0, 0.0}; // outside state bounds

        OsqpMpc controller(sys, cfg);
        auto u = controller.solve(x0);
        REQUIRE(u.has_value());
    }

    SECTION("hard state constraints make out-of-bounds initial state infeasible") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .x_min = (Eigen::Vector2d() << -5.0, -5.0).finished(),
            .x_max = (Eigen::Vector2d() << 5.0, 5.0).finished(),
            .hard_state_constraints = true,
        };

        Eigen::Vector2d x0{10.0, 0.0}; // violates hard state bounds

        OsqpMpc controller(sys, cfg);
        auto result = controller.solve(x0);
        CHECK_FALSE(result.has_value());
    }

    SECTION("warm-starting reduces iterations on second solve") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};

        // First solve (cold)
        auto r1 = controller.solve(x0);
        REQUIRE(r1.has_value());
        int iter1 = controller.diagnostics().iterations;

        // Simulate one step forward so state changes slightly
        Eigen::Vector2d x1 = sys.A * x0 + sys.B * r1.value();

        // Second solve with nearby state (warm-started from first)
        auto r2 = controller.solve(x1);
        REQUIRE(r2.has_value());
        int iter2 = controller.diagnostics().iterations;

        CHECK(iter2 <= iter1);
    }

    SECTION("diagnostics populated after successful solve") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto diag = controller.diagnostics();
        CHECK(diag.iterations > 0);
        CHECK(diag.solve_time > 0.0);
        CHECK(diag.cost >= 0.0);
    }

    SECTION("trajectory returns consistent states and inputs") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto [states, inputs] = controller.trajectory();
        REQUIRE(states.size() == static_cast<std::size_t>(N + 1));
        REQUIRE(inputs.size() == static_cast<std::size_t>(N));

        // states[0] should match x0
        CHECK_THAT(states[0](0), WithinAbs(x0(0), 1e-3));
        CHECK_THAT(states[0](1), WithinAbs(x0(1), 1e-3));

        // States should propagate approximately via dynamics
        for (int k = 0; k < N; ++k) {
            Eigen::Vector2d x_next = sys.A * states[static_cast<std::size_t>(k)]
                + sys.B * inputs[static_cast<std::size_t>(k)];
            CHECK_THAT(states[static_cast<std::size_t>(k + 1)](0),
                       WithinAbs(x_next(0), 1e-2));
            CHECK_THAT(states[static_cast<std::size_t>(k + 1)](1),
                       WithinAbs(x_next(1), 1e-2));
        }
    }

    SECTION("rate constraints limit control change between solves") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << -5.0).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << 5.0).finished(),
            .du_max = (Eigen::Matrix<double, 1, 1>() << 0.2).finished(),
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x{5.0, 0.0};
        double u_prev = 0.0; // initial u_prev_ is zero in mpc

        for (int step = 0; step < 15; ++step) {
            auto u = controller.solve(x);
            REQUIRE(u.has_value());
            double u_cur = (*u)(0);

            CHECK(std::abs(u_cur - u_prev) <= 0.2 + 1e-2);

            u_prev = u_cur;
            x = sys.A * x + sys.B * u.value();
        }
    }

    SECTION("DARE default Qf solves without explicit Qf") {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            // Qf not specified -- DARE should compute it internally
        };

        OsqpMpc controller(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto diag = controller.diagnostics();
        CHECK(diag.status == ctrlpp::solve_status::optimal);
    }
}
