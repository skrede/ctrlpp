#include "bench_metrics.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"
#include "ctrlpp/detail/numerical_diff.h"

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>
#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

// ---------------------------------------------------------------------------
// Dynamics definitions
// ---------------------------------------------------------------------------

constexpr double di2_dt = 0.1;
auto double_integrator_2 = [](const Eigen::Vector2d& x,
                              const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    return Eigen::Vector2d{x(0) + di2_dt * x(1), x(1) + di2_dt * u(0)};
};

auto pendulum_2 = [](const Eigen::Vector2d& x,
                     const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    constexpr double dt = 0.05;
    constexpr double g = 9.81;
    constexpr double l = 1.0;
    double theta = x(0);
    double omega = x(1);
    double alpha = -g / l * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + dt * omega, omega + dt * alpha};
};

constexpr double di4_dt = 0.1;
auto double_integrator_4 = [](const Eigen::Vector4d& x,
                              const Eigen::Vector2d& u) -> Eigen::Vector4d
{
    return Eigen::Vector4d{
        x(0) + di4_dt * x(1),
        x(1) + di4_dt * u(0),
        x(2) + di4_dt * x(3),
        x(3) + di4_dt * u(1)};
};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

using LbfgsbSolver = ctrlpp::argmin_solver<double, ctrlpp::argmin_lbfgsb, false>;
using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, false>;

// ---------------------------------------------------------------------------
// Single-shooting NLP problem builder
// ---------------------------------------------------------------------------

/// Build a single-shooting NLP for box-constrained solvers.
/// Decision variable: u_0, ..., u_{N-1} (N*NU total).
/// Objective: forward-simulate from x0, accumulate quadratic cost.
/// No equality constraints (continuity enforced by forward simulation).
template <std::size_t NX, std::size_t NU, typename Dynamics>
auto build_single_shooting_problem(
    Dynamics dynamics,
    const Eigen::Matrix<double, NX, 1>& x0,
    int horizon,
    const Eigen::Matrix<double, NX, NX>& Q,
    const Eigen::Matrix<double, NU, NU>& R,
    double u_min,
    double u_max) -> ctrlpp::nlp_problem<double>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    int n_vars = horizon * nu;

    // Cost: forward simulate and accumulate quadratic cost
    auto cost_fn = [=](std::span<const double> z) -> double
    {
        double total = 0.0;
        Eigen::Matrix<double, NX, 1> xk = x0;

        for(int k = 0; k < horizon; ++k)
        {
            Eigen::Map<const Eigen::Matrix<double, NU, 1>> uk(z.data() + k * nu);
            total += 0.5 * xk.dot(Q * xk) + 0.5 * uk.dot(R * uk);
            xk = dynamics(xk, uk);
        }

        // Terminal cost
        total += 0.5 * xk.dot(Q * xk);
        return total;
    };

    std::function<double(std::span<const double>)> cost = cost_fn;

    // Gradient via finite differences
    std::function<void(std::span<const double>, std::span<double>)> gradient =
        [cost](std::span<const double> z, std::span<double> grad)
    {
        ctrlpp::detail::finite_diff_gradient<double>(cost, z, grad);
    };

    // No constraints
    std::function<void(std::span<const double>, std::span<double>)> constraints =
        [](std::span<const double>, std::span<double>) {};

    // Variable bounds: input bounds on all decision variables
    Eigen::VectorXd x_lower = Eigen::VectorXd::Constant(n_vars, u_min);
    Eigen::VectorXd x_upper = Eigen::VectorXd::Constant(n_vars, u_max);

    // No constraint bounds (n_constraints = 0)
    Eigen::VectorXd c_lower;
    Eigen::VectorXd c_upper;

    return ctrlpp::nlp_problem<double>{
        .n_vars = n_vars,
        .n_constraints = 0,
        .cost = std::move(cost),
        .gradient = std::move(gradient),
        .constraints = std::move(constraints),
        .x_lower = std::move(x_lower),
        .x_upper = std::move(x_upper),
        .c_lower = std::move(c_lower),
        .c_upper = std::move(c_upper)};
}

// ---------------------------------------------------------------------------
// Benchmark runner
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_benchmark(const std::string& system_name,
                   Dynamics dynamics,
                   int horizon,
                   ankerl::nanobench::Bench& bench,
                   std::ostream& quality_csv)
{
    constexpr int nu = static_cast<int>(NU);

    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;

    auto Q = Eigen::Matrix<double, NX, NX>::Identity();
    auto R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1;

    constexpr double u_min = -10.0;
    constexpr double u_max = 10.0;

    int n_vars = horizon * nu;

    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(horizon);

    // L-BFGS-B solver
    {
        auto problem = build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, u_min, u_max);

        LbfgsbSolver solver{};

        bench.warmup(50).minEpochIterations(50).title(title)
            .run("lbfgsb",
                 [&]
                 {
                     solver.setup(problem);
                     ctrlpp::nlp_update<double> update{};
                     update.x0 = Eigen::VectorXd::Zero(n_vars);
                     auto result = solver.solve(update);
                     ankerl::nanobench::doNotOptimizeAway(result);
                 });

        // Quality: single clean solve
        solver.setup(problem);
        ctrlpp::nlp_update<double> update{};
        update.x0 = Eigen::VectorXd::Zero(n_vars);
        auto result = solver.solve(update);
        auto qm = compute_quality_metrics(problem, result);

        write_quality_csv_row(quality_csv, system_name, "argmin", "lbfgsb",
                              "cold", static_cast<int>(NX), horizon, qm);
    }

    // Argmin SLSQP on same single-shooting problem
    {
        auto problem = build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, u_min, u_max);

        ArgminSlsqp solver{};

        bench.run("argmin_slsqp",
                  [&]
                  {
                      solver.setup(problem);
                      ctrlpp::nlp_update<double> update{};
                      update.x0 = Eigen::VectorXd::Zero(n_vars);
                      auto result = solver.solve(update);
                      ankerl::nanobench::doNotOptimizeAway(result);
                  });

        solver.setup(problem);
        ctrlpp::nlp_update<double> update{};
        update.x0 = Eigen::VectorXd::Zero(n_vars);
        auto result = solver.solve(update);
        auto qm = compute_quality_metrics(problem, result);

        write_quality_csv_row(quality_csv, system_name, "argmin", "slsqp",
                              "cold", static_cast<int>(NX), horizon, qm);
    }

    // NLopt SLSQP on same single-shooting problem
    {
        auto problem = build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, u_min, u_max);

        NloptSolver solver{};

        bench.run("nlopt_slsqp",
                  [&]
                  {
                      solver.setup(problem);
                      ctrlpp::nlp_update<double> update{};
                      update.x0 = Eigen::VectorXd::Zero(n_vars);
                      auto result = solver.solve(update);
                      ankerl::nanobench::doNotOptimizeAway(result);
                  });

        solver.setup(problem);
        ctrlpp::nlp_update<double> update{};
        update.x0 = Eigen::VectorXd::Zero(n_vars);
        auto result = solver.solve(update);
        auto qm = compute_quality_metrics(problem, result);

        write_quality_csv_row(quality_csv, system_name, "nlopt", "slsqp",
                              "cold", static_cast<int>(NX), horizon, qm);
    }
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_lbfgsb_timing.csv");
    std::ofstream quality_csv("bench_lbfgsb_quality.csv");
    write_quality_csv_header(quality_csv);

    // Problem instances (box-constrained single-shooting)
    run_benchmark<2, 1>("double_integrator", double_integrator_2, 10, bench, quality_csv);
    run_benchmark<4, 2>("double_integrator", double_integrator_4, 10, bench, quality_csv);
    run_benchmark<2, 1>("pendulum", pendulum_2, 10, bench, quality_csv);
    run_benchmark<4, 2>("double_integrator", double_integrator_4, 20, bench, quality_csv);

    // Render timing CSV
    bench.render(comma_csv_tpl, timing_csv);
}
