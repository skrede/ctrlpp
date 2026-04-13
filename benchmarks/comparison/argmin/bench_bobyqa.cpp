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
#include <random>
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
// Single-shooting NLP problem builder (box-constrained, no gradient)
// ---------------------------------------------------------------------------

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
    constexpr int nu = static_cast<int>(NU);
    int n_vars = horizon * nu;

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

        total += 0.5 * xk.dot(Q * xk);
        return total;
    };

    std::function<double(std::span<const double>)> cost = cost_fn;

    std::function<void(std::span<const double>, std::span<double>)> gradient =
        [cost](std::span<const double> z, std::span<double> grad)
    {
        ctrlpp::detail::finite_diff_gradient<double>(cost, z, grad);
    };

    std::function<void(std::span<const double>, std::span<double>)> constraints =
        [](std::span<const double>, std::span<double>) {};

    Eigen::VectorXd x_lower = Eigen::VectorXd::Constant(n_vars, u_min);
    Eigen::VectorXd x_upper = Eigen::VectorXd::Constant(n_vars, u_max);

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
// Type aliases
// ---------------------------------------------------------------------------

using BobyqaSolver = ctrlpp::argmin_solver<double, ctrlpp::argmin_bobyqa, false>;
using ArgminCobyla = ctrlpp::argmin_solver<double, ctrlpp::argmin_cobyla>;
using NloptSolver = ctrlpp::nlopt_solver<double>;

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

    // BOBYQA (nablapp, derivative-free, box-constrained)
    {
        auto problem = build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, u_min, u_max);

        ctrlpp::argmin_settings<double> cfg{};
        cfg.max_eval = 1000;
        BobyqaSolver solver{cfg};

        bench.warmup(50).minEpochIterations(20).title(title)
            .run("argmin_bobyqa",
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

        write_quality_csv_row(quality_csv, system_name, "argmin", "bobyqa",
                              "cold", static_cast<int>(NX), horizon, qm);
    }

    // Argmin COBYLA (nablapp, derivative-free, uses constrained bridge with zero constraints)
    {
        auto problem = build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, u_min, u_max);

        ctrlpp::argmin_settings<double> cfg{};
        cfg.max_eval = 1000;
        ArgminCobyla solver{cfg};

        bench.run("argmin_cobyla",
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

        write_quality_csv_row(quality_csv, system_name, "argmin", "cobyla",
                              "cold", static_cast<int>(NX), horizon, qm);
    }

    // NLopt COBYLA (derivative-free baseline)
    {
        auto problem = build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, u_min, u_max);

        ctrlpp::nlopt_settings<double> nlopt_cfg{};
        nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::cobyla;
        nlopt_cfg.max_eval = 1000;
        NloptSolver solver{nlopt_cfg};

        bench.run("nlopt_cobyla",
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

        write_quality_csv_row(quality_csv, system_name, "nlopt", "cobyla",
                              "cold", static_cast<int>(NX), horizon, qm);
    }
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_bobyqa_timing.csv");
    std::ofstream quality_csv("bench_bobyqa_quality.csv");
    write_quality_csv_header(quality_csv);

    // Box-constrained single-shooting problems
    run_benchmark<2, 1>("double_integrator", double_integrator_2, 10, bench, quality_csv);
    run_benchmark<4, 2>("double_integrator", double_integrator_4, 10, bench, quality_csv);
    run_benchmark<2, 1>("pendulum", pendulum_2, 10, bench, quality_csv);
    run_benchmark<4, 2>("double_integrator", double_integrator_4, 20, bench, quality_csv);

    bench.render(comma_csv_tpl, timing_csv);
}
