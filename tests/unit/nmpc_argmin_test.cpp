#include "ctrlpp/mpc/argmin_solver.h"
#ifdef CTRLPP_HAS_NLOPT
#include "ctrlpp/mpc/nlopt_solver.h"
#endif
#include "ctrlpp/nmpc.h"
#ifdef CTRLPP_HAS_OSQP
#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#endif

#include <Eigen/Dense>

#include <utility>

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
        x = double_integrator(x, u->input);
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
        x = double_integrator(x, u->input);
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
        CHECK(u->input(0) >= -0.5 - 1e-6);
        CHECK(u->input(0) <= 0.5 + 1e-6);
        x = double_integrator(x, u->input);
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
        x = pendulum(x, u->input);
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

TEST_CASE("argmin try_setup reports raw-MMA equality rejection as an expected error", "[nmpc][argmin]")
{
    // A minimal NLP with a single equality constraint x0 + x1 = 1. Raw MMA and
    // raw GCMMA cannot represent equality constraints, so try_setup must return
    // the incompatible_equality_constraints error as a value instead of a
    // silently-wrong solve; the auglag-wrapped variant absorbs the equality
    // constraint and must set up successfully on the same problem. Raw MMA is
    // exercised with Constrained=false because the class-body static_assert
    // bars the raw-MMA + constrained-bridge instantiation (the compile-time
    // complement to this runtime reject); the reject still fires by scanning
    // the problem's equalities directly, mirroring nlopt_solver.
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

    SECTION("raw MMA rejects the equality constraint")
    {
        ctrlpp::argmin_solver<double, ctrlpp::argmin_mma, false> solver;
        auto result = solver.try_setup(problem);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::argmin_setup_error::incompatible_equality_constraints);
    }

    SECTION("raw GCMMA rejects the equality constraint")
    {
        ctrlpp::argmin_solver<double, ctrlpp::argmin_gcmma, false> solver;
        auto result = solver.try_setup(problem);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::argmin_setup_error::incompatible_equality_constraints);
    }

    SECTION("auglag-wrapped MMA absorbs the equality constraint")
    {
        ctrlpp::argmin_solver<double, ctrlpp::argmin_auglag<ctrlpp::argmin_mma>> solver;
        auto result = solver.try_setup(problem);
        REQUIRE(result.has_value());
    }
}

TEST_CASE("argmin ftol_rel/xtol_rel drive the relative convergence criteria", "[nmpc][argmin]")
{
    // The settings fields named ftol_rel / xtol_rel must reach argmin's
    // RELATIVE criteria (objective_tolerance_rel / step_tolerance_rel), not the
    // absolute ones. Rosenbrock from (-2, 2) is a sharp lever: a loose relative
    // step tolerance fires step_tolerance_rel_criterion (xtol_reached) within a
    // couple of iterations far short of the optimum, while a tight one runs the
    // full descent to (1, 1). A loose relative objective tolerance likewise
    // terminates objective_tolerance_rel_criterion earlier than the tight run.
    // If the fields were still wired to the absolute setters, neither loosening
    // would change the iterate count (the absolute criteria would stay inert).
    auto make_rosenbrock = []
    {
        ctrlpp::nlp_problem<double> prob;
        prob.n_vars = 2;
        prob.n_constraints = 0;
        prob.cost = [](std::span<const double> x)
        { return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]); };
        prob.gradient = [](std::span<const double> x, std::span<double> g)
        {
            g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
            g[1] = 200.0 * (x[1] - x[0] * x[0]);
        };
        prob.x_lower = Eigen::Vector2d::Constant(-10.0);
        prob.x_upper = Eigen::Vector2d::Constant(10.0);
        prob.c_lower = Eigen::VectorXd{};
        prob.c_upper = Eigen::VectorXd{};
        return prob;
    };

    auto solve_with = [&](double ftol_rel, double xtol_rel) -> ctrlpp::nlp_result<double>
    {
        auto prob = make_rosenbrock();
        ctrlpp::argmin_settings<double> settings{};
        settings.ftol_rel = ftol_rel;
        settings.xtol_rel = xtol_rel;
        settings.max_eval = 500;

        ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false> solver{settings};
        solver.setup(prob);

        ctrlpp::nlp_update<double> update;
        update.x0 = Eigen::Vector2d{-2.0, 2.0};
        return solver.solve(update);
    };

    auto tight = solve_with(1e-12, 1e-12);
    auto loose_xtol = solve_with(1e-12, 5e-1);
    auto loose_ftol = solve_with(5e-1, 1e-12);

    // The tight run must actually reach the Rosenbrock optimum at (1, 1).
    REQUIRE(tight.status == ctrlpp::solve_status::optimal);
    CHECK_THAT(tight.x(0), WithinAbs(1.0, 1e-3));
    CHECK_THAT(tight.x(1), WithinAbs(1.0, 1e-3));

    // A loose relative step tolerance stops early via step_tolerance_rel_criterion.
    CHECK(loose_xtol.status == ctrlpp::solve_status::optimal);
    CHECK(loose_xtol.iterations < tight.iterations);

    // A loose relative objective tolerance stops earlier than the tight descent.
    CHECK(loose_ftol.iterations < tight.iterations);
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
    x = double_integrator(x, u1->input);

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
        x = double_integrator(x, u->input);

        double error = (x - refs[1]).norm();
        max_error = std::max(max_error, error);
    }

    CHECK(max_error < 2.0);
}

// --- MAJ-06 / D-C move & copy safety -----------------------------------------
//
// argmin's solver_core caches the problem BY REFERENCE (const Problem*,
// solver_core.h:657) and the ctrlpp bridge caches &m_problem, so a naive move of
// an nmpc / argmin_solver would relocate the pointed-to bridge and problem while
// the cached back-pointers kept the stale addresses. These cases relocate a
// controller and then drive a full solve through those back-pointers; they are
// compiled into an -fsanitize=address variant (nmpc_argmin_asan_test) so a
// dangling read surfaces as a heap-use-after-free rather than a silent pass.

TEST_CASE("nmpc argmin survives move-then-solve", "[nmpc][argmin][move-safety]")
{
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI source{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};
    const double initial_norm = x.norm();

    // Solve once so the solver variant is emplaced and caches &bridge (layer 1)
    // before the relocation, in addition to the bridge->m_problem link (layer 2).
    auto u0 = source.solve(x);
    REQUIRE(u0.has_value());
    x = double_integrator(x, u0->input);

    // Relocate the controller. A defaulted move must keep both argmin
    // back-pointers valid because the bridge and the problem are heap-stable.
    NmpcDI moved{std::move(source)};

    for(int step = 0; step < 50; ++step)
    {
        auto u = moved.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, u->input);
    }

    REQUIRE(x.norm() < 0.1 * initial_norm);
}

TEST_CASE("nmpc argmin copy is an independent fork", "[nmpc][argmin][move-safety]")
{
    auto config = make_config(10);
    config.Q = 10.0 * Eigen::Matrix2d::Identity();
    config.R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();

    NmpcDI original{double_integrator, config};

    // Warm the original so its formulation state / warm-start is non-trivial;
    // the fork must snapshot, not alias, that state.
    REQUIRE(original.solve(Eigen::Vector2d{1.0, 0.0}).has_value());

    NmpcDI fork{original};

    // Drive both controllers from DIFFERENT initial states in lockstep. If the
    // fork aliased the source's shared state / problem / solver, the interleaved
    // writes would corrupt one another; independence keeps both convergent.
    Eigen::Vector2d xa{1.0, 0.0};
    Eigen::Vector2d xb{-2.0, 1.0};
    const double na0 = xa.norm();
    const double nb0 = xb.norm();

    for(int step = 0; step < 60; ++step)
    {
        auto ua = original.solve(xa);
        auto ub = fork.solve(xb);
        REQUIRE(ua.has_value());
        REQUIRE(ub.has_value());
        xa = double_integrator(xa, ua->input);
        xb = double_integrator(xb, ub->input);
    }

    CHECK(xa.norm() < 0.1 * na0);
    CHECK(xb.norm() < 0.1 * nb0);
}

TEST_CASE("argmin_solver survives move-then-solve", "[argmin][move-safety]")
{
    // The external problem outlives both solvers (stack-local for the whole
    // case), so the only relocation under test is the solver's own bridge.
    ctrlpp::nlp_problem<double> prob;
    prob.n_vars = 2;
    prob.n_constraints = 0;
    prob.cost = [](std::span<const double> x)
    { return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]); };
    prob.gradient = [](std::span<const double> x, std::span<double> g)
    {
        g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
        g[1] = 200.0 * (x[1] - x[0] * x[0]);
    };
    prob.x_lower = Eigen::Vector2d::Constant(-10.0);
    prob.x_upper = Eigen::Vector2d::Constant(10.0);
    prob.c_lower = Eigen::VectorXd{};
    prob.c_upper = Eigen::VectorXd{};

    ctrlpp::argmin_settings<double> settings{};
    settings.max_eval = 500;

    ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false> source{settings};
    source.setup(prob);

    ctrlpp::nlp_update<double> update;
    update.x0 = Eigen::Vector2d{-1.2, 1.0};

    // Solve once so solver_ is emplaced and caches &bridge_ (layer 1) before move.
    auto first = source.solve(update);
    CHECK(first.status != ctrlpp::solve_status::error);

    ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false> moved{std::move(source)};
    auto result = moved.solve(update);

    REQUIRE(result.status == ctrlpp::solve_status::optimal);
    CHECK_THAT(result.x(0), WithinAbs(1.0, 1e-3));
    CHECK_THAT(result.x(1), WithinAbs(1.0, 1e-3));
}

#ifdef CTRLPP_HAS_OSQP
TEST_CASE("mpc osqp survives move-then-solve", "[mpc][osqp][move-safety]")
{
    // Assumption A2: OSQP copies the problem data into its own workspace and
    // caches no pointer into mpc's members (the qp_problem handed to setup is a
    // local in build_initial_qp), so mpc is move-safe. This case pins that.
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt, dt;
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{A, B, C, D};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> source{sys, cfg};
    REQUIRE(source.solve(Eigen::Vector2d{1.0, 0.0}).has_value());

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> moved{std::move(source)};

    Eigen::Vector2d x{1.0, 0.0};
    for(int step = 0; step < 50; ++step)
    {
        auto u = moved.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value();
    }

    CHECK(x.norm() < 0.1);
}
#endif
