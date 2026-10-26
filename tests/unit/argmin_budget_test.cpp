#include "ctrlpp/nmpc.h"

#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <limits>
#include <utility>
#include <variant>
#include <type_traits>

namespace
{

using step_solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false>;
using timed_solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false, Eigen::Dynamic, Eigen::Dynamic,
                                           ctrlpp::argmin_step_and_time_budget>;

template <typename Settings>
concept has_max_time = requires(Settings settings) { settings.max_time; };

template <typename Budget>
concept accepted_budget = requires {
    typename ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false, Eigen::Dynamic, Eigen::Dynamic, Budget>;
};

static_assert(std::is_same_v<ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>,
                             ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, true, Eigen::Dynamic, Eigen::Dynamic,
                                                   ctrlpp::argmin_step_budget>>);

static_assert(std::is_same_v<step_solver::settings_type, ctrlpp::argmin_settings<double>>);
static_assert(!has_max_time<step_solver::settings_type>);
static_assert(std::variant_size_v<step_solver::solver_storage_type> == 2);
static_assert(std::is_same_v<std::variant_alternative_t<1, step_solver::solver_storage_type>,
                             step_solver::step_solver_type>);

static_assert(std::is_same_v<timed_solver::settings_type, ctrlpp::argmin_timed_settings<double>>);
static_assert(std::is_base_of_v<ctrlpp::argmin_settings<double>, timed_solver::settings_type>);
static_assert(has_max_time<timed_solver::settings_type>);
static_assert(std::is_same_v<std::variant_alternative_t<2, timed_solver::solver_storage_type>,
                             timed_solver::timed_solver_type>);

static_assert(accepted_budget<ctrlpp::argmin_step_budget>);
static_assert(accepted_budget<ctrlpp::argmin_step_and_time_budget>);
static_assert(!accepted_budget<int>);

auto make_rosenbrock() -> ctrlpp::nlp_problem<double>
{
    ctrlpp::nlp_problem<double> prob;
    prob.n_vars = 2;
    prob.n_constraints = 0;
    prob.cost = [](std::span<const double> x) {
        return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
    };
    prob.gradient = [](std::span<const double> x, std::span<double> g) {
        g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
        g[1] = 200.0 * (x[1] - x[0] * x[0]);
    };
    prob.x_lower = Eigen::Vector2d::Constant(-10.0);
    prob.x_upper = Eigen::Vector2d::Constant(10.0);
    prob.c_lower = Eigen::VectorXd{};
    prob.c_upper = Eigen::VectorXd{};
    return prob;
}

auto make_constrained_rosenbrock() -> ctrlpp::nlp_problem<double>
{
    auto prob = make_rosenbrock();
    prob.n_constraints = 2;
    prob.constraints = [](std::span<const double> x, std::span<double> c) {
        c[0] = x[0] + x[1];
        c[1] = x[0] - x[1];
    };
    prob.constraint_jacobian = [](std::span<const double>, std::span<double> jac) {
        jac[0] = 1.0;
        jac[1] = 1.0;
        jac[2] = 1.0;
        jac[3] = -1.0;
    };
    prob.c_lower = Eigen::Vector2d::Constant(-std::numeric_limits<double>::infinity());
    prob.c_upper = Eigen::Vector2d::Constant(1.0);
    return prob;
}

struct double_integrator
{
    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        constexpr double dt = 0.1;
        return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)};
    }
};

template <typename Solver>
auto solve_rosenbrock(Solver solver) -> ctrlpp::nlp_result<double>
{
    const auto prob = make_rosenbrock();
    REQUIRE(solver.setup(prob).has_value());
    ctrlpp::nlp_update<double> update;
    update.x0 = Eigen::Vector2d{-2.0, 2.0};
    return solver.solve(update);
}

}

// Only the time-budgeted driver measures the solve, so a zero solve_time names the step-budget driver.
TEST_CASE("argmin step budget runs the step driver", "[argmin][budget]")
{
    const auto result = solve_rosenbrock(step_solver{});

    REQUIRE(result.status == ctrlpp::solve_status::optimal);
    CHECK(result.solve_time == 0.0);
}

TEST_CASE("argmin step-and-time budget selects its driver by max_time", "[argmin][budget]")
{
    const auto stepped = solve_rosenbrock(step_solver{});
    REQUIRE(stepped.status == ctrlpp::solve_status::optimal);

    SECTION("a zero max_time runs the step driver")
    {
        const auto result = solve_rosenbrock(timed_solver{});

        CHECK(result.solve_time == 0.0);
        CHECK(result.iterations == stepped.iterations);
        CHECK(result.x == stepped.x);
    }

    // A minute cannot bind on a two-variable solve, so the timed driver walks the step driver's iterates.
    SECTION("a positive max_time runs the timed driver")
    {
        timed_solver::settings_type settings{};
        settings.max_time = 60.0;

        const auto result = solve_rosenbrock(timed_solver{settings});

        CHECK(result.status == ctrlpp::solve_status::optimal);
        CHECK(result.solve_time > 0.0);
        CHECK(result.iterations == stepped.iterations);
        CHECK(result.x == stepped.x);
    }
}

TEST_CASE("the unconstrained argmin bridge refuses a problem with constraint rows", "[argmin][setup][hardening]")
{
    step_solver solver;
    const auto constrained = make_constrained_rosenbrock();

    const auto refused = solver.setup(constrained);

    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::argmin_setup_error::constraints_unsupported);

    SECTION("a problem without constraint rows still sets up and solves")
    {
        const auto result = solve_rosenbrock(std::move(solver));

        REQUIRE(result.status == ctrlpp::solve_status::optimal);
        CHECK_THAT(result.x(0), Catch::Matchers::WithinAbs(1.0, 1e-3));
        CHECK_THAT(result.x(1), Catch::Matchers::WithinAbs(1.0, 1e-3));
    }
}

// The transcription always poses the dynamics as equality rows, so no NMPC problem fits the unconstrained bridge.
TEST_CASE("nmpc on the unconstrained argmin bridge is refused at construction", "[argmin][setup][hardening]")
{
    using controller_type = ctrlpp::nmpc_dynamic<double, 2, 1, step_solver, double_integrator>;
    const ctrlpp::nmpc_config<double, 2, 1> config{
        .horizon = 5,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };

    const auto controller = controller_type::create(double_integrator{}, config);

    REQUIRE_FALSE(controller.has_value());
    const auto* cause = std::get_if<ctrlpp::argmin_setup_error>(&controller.error());
    REQUIRE(cause != nullptr);
    CHECK(*cause == ctrlpp::argmin_setup_error::constraints_unsupported);
}
