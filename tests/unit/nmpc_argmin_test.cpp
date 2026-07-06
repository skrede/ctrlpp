#include "ctrlpp/mpc/argmin_solver.h"
#ifdef CTRLPP_HAS_NLOPT
#include "ctrlpp/mpc/nlopt_solver.h"
#endif
#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <cmath>
#include <vector>
#include <cstddef>
#include <type_traits>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d { return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

auto pendulum = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    constexpr double pdt = 0.05;
    constexpr double g = 9.81;
    constexpr double l = 1.0;
    double theta = x(0);
    double omega = x(1);
    double alpha = -g / l * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + pdt * omega, omega + pdt * alpha};
};

auto make_config(int horizon = 10) -> ctrlpp::nmpc_config<double, NX, NU>
{
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
    };
}

using ArgminSolver = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;
using NmpcDI = ctrlpp::nmpc<double, NX, NU, ArgminSolver, decltype(double_integrator)>;
using NmpcPend = ctrlpp::nmpc<double, NX, NU, ArgminSolver, decltype(pendulum)>;

}

TEST_CASE("argmin_solver satisfies nlp_solver concept", "[argmin]")
{
    static_assert(ctrlpp::nlp_solver<ArgminSolver>, "argmin_solver<double, argmin_slsqp> must satisfy nlp_solver concept");
    static_assert(std::is_same_v<decltype(ctrlpp::nlp_result<double>{}.status), ctrlpp::solve_status>, "nlp_result::status must stay a plain solve_status value");
}

TEST_CASE("argmin_solver satisfies nlp_stepper concept", "[argmin]")
{
    static_assert(ctrlpp::nlp_stepper<ArgminSolver>, "argmin_solver<double, argmin_slsqp> must satisfy nlp_stepper concept");
}

TEST_CASE("nmpc argmin regulation", "[nmpc][argmin]")
{
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    double initial_norm = x.norm();

    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    REQUIRE(x.norm() < 0.1 * initial_norm);
}

TEST_CASE("nmpc argmin setpoint tracking", "[nmpc][argmin]")
{
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.0, 0.0};
    Eigen::Vector2d x_ref{2.0, 0.0};

    for(int step = 0; step < 80; ++step)
    {
        auto u = controller.solve(x, x_ref);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    REQUIRE((x - x_ref).norm() < 0.5);
}

TEST_CASE("nmpc argmin input box constraints", "[nmpc][argmin]")
{
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();
    config.u_min = Eigen::Matrix<double, 1, 1>{-0.5};
    config.u_max = Eigen::Matrix<double, 1, 1>{0.5};

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{5.0, 0.0};

    for(int step = 0; step < 30; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        CHECK((*u)(0) >= -0.5 - 1e-6);
        CHECK((*u)(0) <= 0.5 + 1e-6);
        x = double_integrator(x, *u);
    }
}

TEST_CASE("nmpc argmin pendulum regulation", "[nmpc][argmin]")
{
    auto config = make_config(5);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcPend controller{pendulum, config};

    Eigen::Vector2d x{0.3, 0.0};

    for(int step = 0; step < 20; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = pendulum(x, *u);
    }

    REQUIRE(x.norm() < 0.3);
}

TEST_CASE("nmpc argmin all policies compile", "[nmpc][argmin]")
{
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_slsqp>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_nw_sqp>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_auglag<>>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_cobyla>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_isres>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_lbfgsb, false>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_byrd_lbfgsb, false>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_bobyqa, false>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_mma, false>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_gcmma, false>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_auglag<ctrlpp::argmin_mma>>>);
    static_assert(ctrlpp::nlp_solver<ctrlpp::argmin_solver<double, ctrlpp::argmin_auglag<ctrlpp::argmin_gcmma>>>);
}

#ifdef CTRLPP_HAS_NLOPT
TEST_CASE("nlopt auglag_eq + ld_mma smoke", "[nmpc][argmin][nlopt]")
{
    // Drive the NLopt AUGLAG_EQ + LD_MMA composition directly on the
    // NMPC-built nlp_problem. The smoke gate is solver-level: setup() must
    // wire the inner local optimizer before optimize() is called, no
    // exceptions are thrown, the solve returns a non-error status, and
    // the solution vector is dimensionally correct with finite entries.
    //
    // We deliberately bypass the nmpc<>::solve optional-collapsing layer
    // here: AUGLAG_EQ shares the outer max_eval budget with the inner
    // LD_MMA local-solve (per NLopt semantics, evaluations are counted
    // jointly), so a tight budget routinely terminates with status
    // max_iterations rather than optimal. The bench layer documents this
    // ambiguity in its CSV header; the unit-level gate here is just that
    // the composition runs end-to-end without throwing or erroring.

    using NloptSolver = ctrlpp::nlopt_solver<double>;
    using NmpcDIN = ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator)>;

    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDIN controller{double_integrator, config};
    const auto& problem = controller.problem();

    ctrlpp::nlopt_settings<double> nlopt_cfg{};
    nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::auglag_mma;
    nlopt_cfg.ftol_rel = 1e-6;
    nlopt_cfg.xtol_rel = 1e-6;
    nlopt_cfg.max_eval = 500;
    nlopt_cfg.constraint_tol = 1e-6;

    NloptSolver solver{nlopt_cfg};
    REQUIRE_NOTHROW(solver.setup(problem));

    ctrlpp::nlp_update<double> update;
    update.x0 = Eigen::VectorXd::Zero(problem.n_vars);

    ctrlpp::nlp_result<double> result;
    REQUIRE_NOTHROW(result = solver.solve(update));

    CHECK(result.status != ctrlpp::solve_status::error);
    CHECK(result.x.size() == problem.n_vars);
    CHECK(std::isfinite(result.objective));
    for(Eigen::Index i = 0; i < result.x.size(); ++i)
        CHECK(std::isfinite(result.x[i]));
}

TEST_CASE("nlopt try_setup reports equality-constraint rejection as an expected error", "[nmpc][argmin][nlopt]")
{
    // A minimal NLP with a single equality constraint x0 + x1 = 1. Raw MMA and
    // raw CCSAQ cannot handle equality constraints, so try_setup must return
    // the incompatible_equality_constraints error as a value instead of
    // throwing; the auglag-wrapped variant absorbs the equality constraint and
    // must set up successfully on the same problem.
    ctrlpp::nlp_problem<double> problem{};
    problem.n_vars = 2;
    problem.n_constraints = 1;
    problem.cost = [](std::span<const double> x) { return x[0] * x[0] + x[1] * x[1]; };
    problem.gradient = [](std::span<const double> x, std::span<double> g)
    {
        g[0] = 2.0 * x[0];
        g[1] = 2.0 * x[1];
    };
    problem.constraints = [](std::span<const double> x, std::span<double> c) { c[0] = x[0] + x[1]; };
    problem.c_lower = Eigen::VectorXd::Constant(1, 1.0);
    problem.c_upper = Eigen::VectorXd::Constant(1, 1.0);

    SECTION("mma rejects the equality constraint")
    {
        ctrlpp::nlopt_settings<double> settings{};
        settings.algorithm = ctrlpp::nlopt_algorithm::mma;

        ctrlpp::nlopt_solver<double> solver{settings};
        auto result = solver.try_setup(problem);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::nlopt_setup_error::incompatible_equality_constraints);
    }

    SECTION("ccsaq rejects the equality constraint")
    {
        ctrlpp::nlopt_settings<double> settings{};
        settings.algorithm = ctrlpp::nlopt_algorithm::ccsaq;

        ctrlpp::nlopt_solver<double> solver{settings};
        auto result = solver.try_setup(problem);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::nlopt_setup_error::incompatible_equality_constraints);
    }

    SECTION("auglag_mma absorbs the equality constraint")
    {
        ctrlpp::nlopt_settings<double> settings{};
        settings.algorithm = ctrlpp::nlopt_algorithm::auglag_mma;

        ctrlpp::nlopt_solver<double> solver{settings};
        auto result = solver.try_setup(problem);
        REQUIRE(result.has_value());
    }
}
#endif

TEST_CASE("nmpc argmin warm-start benefit", "[nmpc][argmin]")
{
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

TEST_CASE("nmpc argmin trajectory tracking", "[nmpc][argmin]")
{
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI controller{double_integrator, config};

    Eigen::Vector2d x{0.0, 0.0};
    double max_error = 0.0;

    for(int step = 0; step < 50; ++step)
    {
        double t = step * dt;

        std::vector<Eigen::Vector2d> refs;
        refs.reserve(11);
        for(int k = 0; k <= 10; ++k)
        {
            double tk = t + k * dt;
            refs.push_back(Eigen::Vector2d{std::sin(0.5 * tk), 0.5 * std::cos(0.5 * tk)});
        }

        auto u = controller.solve(x, std::span<const Eigen::Vector2d>{refs});
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);

        double error = (x - refs[1]).norm();
        max_error = std::max(max_error, error);
    }

    CHECK(max_error < 2.0);
}
