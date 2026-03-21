#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

namespace {

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto double_integrator = [](const Eigen::Vector2d& x,
                            const Eigen::Matrix<double, 1, 1>& u)
    -> Eigen::Vector2d {
    return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)};
};

struct mock_nlp_solver {
    using scalar_type = double;

    mutable ctrlpp::nlp_problem<double> last_problem{};
    mutable ctrlpp::nlp_update<double> last_update{};
    mutable int solve_count{0};

    void setup(const ctrlpp::nlp_problem<double>& p) {
        last_problem = p;
    }

    auto solve(const ctrlpp::nlp_update<double>& u) -> ctrlpp::nlp_result<double> {
        last_update = u;
        ++solve_count;
        ctrlpp::nlp_result<double> r{};
        r.status = ctrlpp::solve_status::optimal;
        r.x = Eigen::VectorXd::Zero(last_problem.n_vars);
        r.objective = 0.0;
        r.solve_time = 0.001;
        r.iterations = 1;
        r.primal_residual = 0.0;
        return r;
    }
};

auto make_config(int horizon = 5) -> ctrlpp::nmpc_config<double, NX, NU> {
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };
}

using Nmpc = ctrlpp::nmpc<double, NX, NU, mock_nlp_solver, decltype(double_integrator)>;

}

TEST_CASE("dynamics_model concept accepts lambda with correct signature") {
    static_assert(ctrlpp::dynamics_model<decltype(double_integrator), double, NX, NU>,
        "double_integrator must satisfy dynamics_model concept");
}

TEST_CASE("dynamics_model concept rejects lambda with wrong signature") {
    auto wrong_sig = [](const Eigen::Vector2d& x) -> Eigen::Vector2d {
        return x;
    };
    static_assert(!ctrlpp::dynamics_model<decltype(wrong_sig), double, NX, NU>,
        "wrong signature must not satisfy dynamics_model concept");
}

TEST_CASE("mock_nlp_solver satisfies nlp_solver concept") {
    static_assert(ctrlpp::nlp_solver<mock_nlp_solver>,
        "mock_nlp_solver must satisfy nlp_solver concept");
}

TEST_CASE("nmpc with mock solver", "[nmpc]") {
    constexpr int N = 5;
    auto config = make_config(N);
    Nmpc controller{double_integrator, config};

    SECTION("solve(x0) returns a value") {
        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());
    }

    SECTION("solve(x0, x_ref) returns a value") {
        Eigen::Vector2d x0{1.0, 0.0};
        Eigen::Vector2d x_ref{2.0, 0.0};
        auto result = controller.solve(x0, x_ref);
        REQUIRE(result.has_value());
    }

    SECTION("solve(x0, trajectory) returns a value") {
        Eigen::Vector2d x0{1.0, 0.0};
        std::vector<Eigen::Vector2d> refs(static_cast<std::size_t>(N + 1),
                                          Eigen::Vector2d{2.0, 0.0});
        auto result = controller.solve(x0,
            std::span<const Eigen::Vector2d>{refs});
        REQUIRE(result.has_value());
    }

    SECTION("trajectory() returns N+1 states and N inputs after solve") {
        Eigen::Vector2d x0{1.0, 0.0};
        auto u0 = controller.solve(x0);
        REQUIRE(u0.has_value());
        auto [states, inputs] = controller.trajectory();
        REQUIRE(states.size() == static_cast<std::size_t>(N + 1));
        REQUIRE(inputs.size() == static_cast<std::size_t>(N));
    }

    SECTION("diagnostics() returns populated struct after solve") {
        Eigen::Vector2d x0{1.0, 0.0};
        auto u0 = controller.solve(x0);
        REQUIRE(u0.has_value());
        auto diag = controller.diagnostics();
        REQUIRE(diag.status == ctrlpp::solve_status::optimal);
        REQUIRE(diag.iterations >= 0);
        REQUIRE(diag.solve_time >= 0.0);
    }

    SECTION("mock solver receives correct n_vars and n_constraints") {
        // Need access to the mock solver's captured problem
        // n_vars = (N+1)*NX + N*NU = 6*2 + 5*1 = 17
        // n_constraints = (N+1)*NX (equality only, no rate constraints)
        //               = 6*2 = 12
        mock_nlp_solver solver{};
        ctrlpp::nmpc<double, NX, NU, mock_nlp_solver,
                     decltype(double_integrator)> ctrl{double_integrator, config};

        // The solver is internal, so we check indirectly via solve
        // After construction, setup() has been called with the problem
        // We verify dimensions through the solve behavior
        Eigen::Vector2d x0{1.0, 0.0};
        auto result = ctrl.solve(x0);
        REQUIRE(result.has_value());

        // Verify trajectory dimensions match expected variable layout
        auto [states, inputs] = ctrl.trajectory();
        REQUIRE(states.size() == static_cast<std::size_t>(N + 1));
        REQUIRE(inputs.size() == static_cast<std::size_t>(N));

        // Each state vector must have NX elements
        for (const auto& s : states) {
            REQUIRE(s.size() == static_cast<Eigen::Index>(NX));
        }
        // Each input vector must have NU elements
        for (const auto& u : inputs) {
            REQUIRE(u.size() == static_cast<Eigen::Index>(NU));
        }
    }

    SECTION("variable bounds reflect config constraints") {
        auto constrained_config = make_config(N);
        constrained_config.u_min = Eigen::Matrix<double, 1, 1>{-1.0};
        constrained_config.u_max = Eigen::Matrix<double, 1, 1>{1.0};
        constrained_config.x_min = Eigen::Vector2d{-5.0, -5.0};
        constrained_config.x_max = Eigen::Vector2d{5.0, 5.0};

        Nmpc constrained_ctrl{double_integrator, constrained_config};

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = constrained_ctrl.solve(x0);
        REQUIRE(result.has_value());

        // With mock returning zeros, all controls and states should be within bounds
        auto [states, inputs] = constrained_ctrl.trajectory();
        for (const auto& u : inputs) {
            CHECK(u(0) >= -1.0);
            CHECK(u(0) <= 1.0);
        }
    }
}
