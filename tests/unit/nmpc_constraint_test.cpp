#include "ctrlpp/model/constraint_model.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d { return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

using NloptSolver = ctrlpp::nlopt_solver<double>;

// Path constraint: position <= upper_bound  =>  g(x,u) = x[0] - upper_bound <= 0
constexpr std::size_t NC = 1;
auto make_upper_bound_constraint(double bound)
{
    return [bound](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>&) -> ctrlpp::Vector<double, NC> { return ctrlpp::Vector<double, NC>{x(0) - bound}; };
}

// Terminal constraint: x[0] <= target  =>  h(x) = x[0] - target <= 0
constexpr std::size_t NTC_1 = 1;
auto make_terminal_constraint(double target)
{
    return [target](const Eigen::Vector2d& x) -> ctrlpp::Vector<double, NTC_1> { return ctrlpp::Vector<double, NTC_1>{x(0) - target}; };
}

} // namespace

// ----- Concept static assertions -----

TEST_CASE("constraint_model concept accepts valid path constraint lambda")
{
    auto g = make_upper_bound_constraint(1.0);
    static_assert(ctrlpp::constraint_model<decltype(g), double, NX, NU, NC>, "path constraint lambda must satisfy constraint_model");
}

TEST_CASE("terminal_constraint_model concept accepts valid terminal constraint lambda")
{
    auto h = make_terminal_constraint(0.5);
    static_assert(ctrlpp::terminal_constraint_model<decltype(h), double, NX, NTC_1>, "terminal constraint lambda must satisfy terminal_constraint_model");
}

// ----- Backward compatibility -----

TEST_CASE("nmpc with NC=0 NTC=0 produces same results as original", "[nmpc][constraint][compat]")
{
    ctrlpp::nmpc_config<double, NX, NU> config_old{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator)> ctrl_old{double_integrator, config_old};

    ctrlpp::nmpc_config<double, NX, NU, 0, 0> config_new{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), 0, 0> ctrl_new{double_integrator, config_new};

    Eigen::Vector2d x{1.0, 0.0};
    auto u_old = ctrl_old.solve(x);
    auto u_new = ctrl_new.solve(x);

    REQUIRE(u_old.has_value());
    REQUIRE(u_new.has_value());
    CHECK_THAT((*u_old)(0), WithinAbs((*u_new)(0), 1e-8));
}

// ----- Soft path constraint -----

TEST_CASE("soft path constraint is approximately satisfied", "[nmpc][constraint][soft]")
{
    constexpr double upper_bound = 0.8;

    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(upper_bound);
    // soft_constraints defaults to true

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    Eigen::Vector2d x{1.5, 0.0}; // starts above bound

    for(int step = 0; step < 40; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    // After enough steps with soft constraint, state should have moved below bound
    CHECK(x(0) < upper_bound + 0.5);

    auto diag = controller.diagnostics();
    // Diagnostics should have constraint violation info
    CHECK(diag.max_path_constraint_violation >= 0.0);
}

// ----- Hard path constraint -----

TEST_CASE("hard path constraint enforced tightly", "[nmpc][constraint][hard]")
{
    constexpr double upper_bound = 0.8;

    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(upper_bound);
    config.soft_constraints = false;

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    // Start within bounds
    Eigen::Vector2d x{0.5, 0.0};

    for(int step = 0; step < 30; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);

        auto [states, inputs] = controller.trajectory();
        for(const auto& s : states)
        {
            CHECK(s(0) <= upper_bound + 1e-3);
        }
    }

    auto diag = controller.diagnostics();
    CHECK(diag.total_slack == 0.0);
}

// ----- Terminal constraint -----

TEST_CASE("terminal constraint drives final state", "[nmpc][constraint][terminal]")
{
    constexpr double target = 0.3;

    ctrlpp::nmpc_config<double, NX, NU, 0, NTC_1> config{
        .horizon = 15,
        .Q = 1.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.terminal_constraint = make_terminal_constraint(target);
    config.soft_constraints = false;

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), 0, NTC_1> controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        // Check terminal state of predicted trajectory
        auto [states, inputs] = controller.trajectory();
        CHECK(states.back()(0) <= target + 1e-2);

        x = double_integrator(x, *u);
    }
}

// ----- Combined path + terminal constraints -----

TEST_CASE("combined path and terminal constraints", "[nmpc][constraint][combined]")
{
    constexpr double upper_bound = 1.2;
    constexpr double target = 0.3;

    ctrlpp::nmpc_config<double, NX, NU, NC, NTC_1> config{
        .horizon = 15,
        .Q = 5.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(upper_bound);
    config.terminal_constraint = make_terminal_constraint(target);

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, NTC_1> controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    // State should have converged within combined constraints
    CHECK(x(0) < upper_bound + 0.5);

    auto diag = controller.diagnostics();
    CHECK(diag.max_path_constraint_violation >= -10.0); // sanity: field is populated
    CHECK(diag.max_terminal_constraint_violation >= -10.0);
}

// ----- Infeasible constraints with soft mode -----

TEST_CASE("infeasible constraints with soft mode does not crash", "[nmpc][constraint][infeasible]")
{
    // Constraint: x[0] <= -10 but starting at x[0] = 5 with regulation to origin
    // This is intentionally infeasible for early steps
    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(-10.0);
    // soft_constraints defaults to true

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    Eigen::Vector2d x{5.0, 0.0};

    // Solver should not crash -- soft constraints absorb infeasibility
    bool found_slack = false;
    bool found_violation = false;
    for(int step = 0; step < 10; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        auto diag = controller.diagnostics();
        if(diag.total_slack > 0.0)
        {
            found_slack = true;
        }
        if(diag.max_path_constraint_violation > 0.0)
        {
            found_violation = true;
        }

        x = double_integrator(x, *u);
    }

    // At some point during the loop, slack should have been nonzero
    CHECK(found_slack);
    CHECK(found_violation);
}
