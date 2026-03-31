#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/nlp_formulation.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{ return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

using NloptSolver = ctrlpp::nlopt_solver<double>;
using NmpcDI = ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator)>;

}

// ---- nlp_solver.h: concept satisfaction static asserts ----

TEST_CASE("nlopt_solver satisfies nlp_solver concept", "[nmpc][coverage][concept]")
{
    static_assert(ctrlpp::nlp_solver<NloptSolver>);
}

TEST_CASE("nlp_solver concept rejects type without setup", "[nmpc][coverage][concept]")
{
    struct no_setup
    {
        using scalar_type = double;
        auto solve(const ctrlpp::nlp_update<double>&) -> ctrlpp::nlp_result<double> { return {}; }
    };
    static_assert(!ctrlpp::nlp_solver<no_setup>);
}

TEST_CASE("nlp_solver concept rejects type without scalar_type", "[nmpc][coverage][concept]")
{
    struct no_scalar
    {
        void setup(const ctrlpp::nlp_problem<double>&) {}
        auto solve(const ctrlpp::nlp_update<double>&) -> ctrlpp::nlp_result<double> { return {}; }
    };
    static_assert(!ctrlpp::nlp_solver<no_scalar>);
}

// ---- nmpc_config.h: various config field combinations ----

TEST_CASE("nmpc_config default_penalty returns 1e4 for non-zero N", "[nmpc][coverage][config]")
{
    auto p = ctrlpp::detail::default_penalty<double, 2>();
    CHECK_THAT(p(0), WithinAbs(1e4, 1e-6));
    CHECK_THAT(p(1), WithinAbs(1e4, 1e-6));
}

TEST_CASE("nmpc_config default_penalty returns empty vector for N=0", "[nmpc][coverage][config]")
{
    auto p = ctrlpp::detail::default_penalty<double, 0>();
    CHECK(p.size() == 0);
}

TEST_CASE("nmpc_config with explicit Qf terminal cost", "[nmpc][coverage][config]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
        .Qf = 5.0 * Eigen::Matrix2d::Identity(),
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());

    auto diag = controller.diagnostics();
    CHECK(diag.status == ctrlpp::solve_status::optimal);
}

TEST_CASE("nmpc_config with all optional bounds set", "[nmpc][coverage][config]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 8,
        .Q = 5.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
        .u_min = Eigen::Matrix<double, 1, 1>{-2.0},
        .u_max = Eigen::Matrix<double, 1, 1>{2.0},
        .x_min = Eigen::Vector2d{-5.0, -5.0},
        .x_max = Eigen::Vector2d{5.0, 5.0},
        .du_max = Eigen::Matrix<double, 1, 1>{0.5},
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{2.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());
    CHECK((*u)(0) >= -2.0 - 1e-4);
    CHECK((*u)(0) <= 2.0 + 1e-4);
}

// ---- nlp_formulation.h: custom stage and terminal cost functions ----

TEST_CASE("custom stage cost with asymmetric weighting", "[nmpc][coverage][nlp]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };

    // Custom: penalize position heavily, ignore velocity entirely
    config.stage_cost = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> double
    { return 50.0 * x(0) * x(0) + 0.01 * u(0) * u(0); };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 1.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());

    // Should produce aggressive control to reduce position
    CHECK((*u)(0) < -0.01);
}

TEST_CASE("custom terminal cost only (stage cost from Q/R)", "[nmpc][coverage][nlp]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    config.terminal_cost = [](const Eigen::Vector2d& x) -> double
    { return 200.0 * x.squaredNorm(); };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());
    CHECK(std::isfinite((*u)(0)));
}

// ---- nlp_formulation.h: rate constraints (du_max) in NLP ----

TEST_CASE("NLP rate constraints limit control change", "[nmpc][coverage][nlp]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
        .du_max = Eigen::Matrix<double, 1, 1>{0.15},
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{3.0, 0.0};
    double u_prev = 0.0;

    for(int step = 0; step < 15; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        double du = std::abs((*u)(0) - u_prev);
        CHECK(du <= 0.15 + 1e-3);

        u_prev = (*u)(0);
        x = double_integrator(x, *u);
    }
}

// ---- nlp_formulation.h: path constraints with slack ----

TEST_CASE("soft path constraint with custom penalty weight", "[nmpc][coverage][nlp]")
{
    constexpr std::size_t NC = 1;

    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    config.path_constraint = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>&)
        -> ctrlpp::Vector<double, NC>
    { return ctrlpp::Vector<double, NC>{x(0) - 0.5}; };

    config.path_penalty = ctrlpp::Vector<double, NC>{1e6};

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    Eigen::Vector2d x{2.0, 0.0};
    for(int step = 0; step < 30; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    // High penalty should drive state toward constraint, verify it converges
    CHECK(x(0) < 5.0);
}

// ---- nmpc.h: warm-start shifting over consecutive solves ----

TEST_CASE("warm-start shift reduces solve time on consecutive calls", "[nmpc][coverage]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    // Cold solve
    auto u1 = controller.solve(x);
    REQUIRE(u1.has_value());
    auto diag1 = controller.diagnostics();

    // Step forward
    x = double_integrator(x, *u1);

    // Warm solve
    auto u2 = controller.solve(x);
    REQUIRE(u2.has_value());
    auto diag2 = controller.diagnostics();

    // Step forward
    x = double_integrator(x, *u2);

    // Third warm solve
    auto u3 = controller.solve(x);
    REQUIRE(u3.has_value());
    auto diag3 = controller.diagnostics();

    // Warm-started solves should use fewer or equal iterations
    CHECK(diag3.iterations <= diag1.iterations);
}

// ---- nmpc.h: trajectory extraction ----

TEST_CASE("trajectory extraction returns dynamically consistent states", "[nmpc][coverage]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 8,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x0{1.5, -0.5};
    auto u = controller.solve(x0);
    REQUIRE(u.has_value());

    auto [states, inputs] = controller.trajectory();
    REQUIRE(states.size() == 9);
    REQUIRE(inputs.size() == 8);

    // Verify dynamics consistency along the trajectory
    for(std::size_t k = 0; k < inputs.size(); ++k)
    {
        Eigen::Vector2d x_pred = double_integrator(states[k], inputs[k]);
        CHECK_THAT(states[k + 1](0), WithinAbs(x_pred(0), 0.1));
        CHECK_THAT(states[k + 1](1), WithinAbs(x_pred(1), 0.1));
    }
}

// ---- nmpc.h: constraint violation diagnostics ----

TEST_CASE("diagnostics report zero constraint violation when unconstrained", "[nmpc][coverage]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 5,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.5, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());

    auto diag = controller.diagnostics();
    CHECK(diag.max_path_constraint_violation == 0.0);
    CHECK(diag.max_terminal_constraint_violation == 0.0);
    CHECK(diag.total_slack == 0.0);
}

// ---- nmpc.h: solve with partial reference trajectory ----

TEST_CASE("solve with shorter-than-horizon reference trajectory", "[nmpc][coverage]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.0, 0.0};

    // Provide only 5 reference points for a horizon of 10
    // The remaining should be filled with the last reference
    std::vector<Eigen::Vector2d> refs(5, Eigen::Vector2d{2.0, 0.0});

    auto u = controller.solve(x, std::span<const Eigen::Vector2d>{refs});
    REQUIRE(u.has_value());

    // Control should push toward the reference
    CHECK((*u)(0) > 0.0);
}

TEST_CASE("solve with exact-length reference trajectory", "[nmpc][coverage]")
{
    constexpr int N = 8;
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = N,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.0, 0.0};
    std::vector<Eigen::Vector2d> refs(static_cast<std::size_t>(N + 1), Eigen::Vector2d{1.0, 0.0});

    auto u = controller.solve(x, std::span<const Eigen::Vector2d>{refs});
    REQUIRE(u.has_value());
    CHECK((*u)(0) > 0.0);
}

// ---- nlp_formulation.h: Qf fallback to Q when not set ----

TEST_CASE("NLP formulation uses Q as terminal cost when Qf is not set", "[nmpc][coverage][nlp]")
{
    ctrlpp::nmpc_config<double, NX, NU> config_no_qf{
        .horizon = 5,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc_config<double, NX, NU> config_with_qf{
        .horizon = 5,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
        .Qf = 10.0 * Eigen::Matrix2d::Identity(), // same as Q
    };

    NmpcDI ctrl_no_qf{double_integrator, config_no_qf};
    NmpcDI ctrl_with_qf{double_integrator, config_with_qf};

    Eigen::Vector2d x{1.0, 0.0};
    auto u1 = ctrl_no_qf.solve(x);
    auto u2 = ctrl_with_qf.solve(x);

    REQUIRE(u1.has_value());
    REQUIRE(u2.has_value());

    // When Qf == Q, both should produce the same control
    CHECK_THAT((*u1)(0), WithinAbs((*u2)(0), 1e-3));
}

// ---- nlp_formulation.h: state bounds in NLP variable bounds ----

TEST_CASE("NLP state bounds are respected in trajectory", "[nmpc][coverage][nlp]")
{
    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
        .x_min = Eigen::Vector2d{-1.5, -1.5},
        .x_max = Eigen::Vector2d{1.5, 1.5},
    };

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.5};

    for(int step = 0; step < 20; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);

        auto [states, inputs] = controller.trajectory();
        for(const auto& s : states)
        {
            CHECK(s(0) >= -1.5 - 1e-3);
            CHECK(s(0) <= 1.5 + 1e-3);
        }
    }
}

// ---- nmpc.h: infeasible solver returns nullopt ----

TEST_CASE("nmpc returns nullopt on solver failure via mock", "[nmpc][coverage]")
{
    struct failing_nlp_solver
    {
        using scalar_type = double;

        mutable ctrlpp::nlp_problem<double> prob{};

        void setup(const ctrlpp::nlp_problem<double>& p) { prob = p; }

        auto solve(const ctrlpp::nlp_update<double>&) -> ctrlpp::nlp_result<double>
        {
            return {
                .status = ctrlpp::solve_status::error,
                .x = Eigen::VectorXd::Zero(prob.n_vars),
                .objective = 0.0,
                .solve_time = 0.0,
                .iterations = 0,
                .primal_residual = 0.0,
            };
        }
    };

    static_assert(ctrlpp::nlp_solver<failing_nlp_solver>);

    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 5,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc<double, NX, NU, failing_nlp_solver, decltype(double_integrator)> controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    CHECK_FALSE(u.has_value());
}

// ---- nmpc.h: solved_inaccurate is accepted ----

TEST_CASE("nmpc accepts solved_inaccurate status", "[nmpc][coverage]")
{
    struct inaccurate_nlp_solver
    {
        using scalar_type = double;

        mutable ctrlpp::nlp_problem<double> prob{};

        void setup(const ctrlpp::nlp_problem<double>& p) { prob = p; }

        auto solve(const ctrlpp::nlp_update<double>&) -> ctrlpp::nlp_result<double>
        {
            return {
                .status = ctrlpp::solve_status::solved_inaccurate,
                .x = Eigen::VectorXd::Zero(prob.n_vars),
                .objective = 0.5,
                .solve_time = 0.001,
                .iterations = 10,
                .primal_residual = 1e-4,
            };
        }
    };

    static_assert(ctrlpp::nlp_solver<inaccurate_nlp_solver>);

    ctrlpp::nmpc_config<double, NX, NU> config{
        .horizon = 5,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc<double, NX, NU, inaccurate_nlp_solver, decltype(double_integrator)> controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    auto u = controller.solve(x);
    REQUIRE(u.has_value());
}
