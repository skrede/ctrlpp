#include "bench_metrics.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>
#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <random>
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
// NMPC config factory
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU>
auto make_nmpc_config(int horizon) -> ctrlpp::nmpc_config<double, NX, NU>
{
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix<double, NX, NX>::Identity(),
        .R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1,
    };
}

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminIsres = ctrlpp::argmin_solver<double, ctrlpp::argmin_isres>;

auto warm_start_label(ctrlpp::warm_start_mode ws) -> std::string
{
    switch(ws)
    {
    case ctrlpp::warm_start_mode::cold:        return "cold";
    case ctrlpp::warm_start_mode::primal_only: return "primal_only";
    case ctrlpp::warm_start_mode::curvature:   return "curvature";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// Benchmark runner
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_benchmark(const std::string& system_name,
                   Dynamics dynamics,
                   int horizon,
                   ankerl::nanobench::Bench& bench,
                   std::ostream& quality_csv,
                   ctrlpp::warm_start_mode ws_mode)
{
    auto config = make_nmpc_config<NX, NU>(horizon);
    auto ws_label = warm_start_label(ws_mode);

    // NLopt solver (always cold -- NLopt has no warm-start)
    ctrlpp::nlopt_settings<double> nlopt_cfg{};
    nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::isres;
    nlopt_cfg.max_eval = 2000;
    NloptSolver nlopt_solver{nlopt_cfg};

    // Argmin solver with requested warm-start
    ctrlpp::argmin_settings<double> argmin_cfg{};
    argmin_cfg.warm_start = ws_mode;
    argmin_cfg.max_eval = 2000;
    ArgminIsres argmin_solver{argmin_cfg};

    ctrlpp::nmpc<double, NX, NU, NloptSolver, Dynamics> nmpc_nlopt{dynamics, config};
    ctrlpp::nmpc<double, NX, NU, ArgminIsres, Dynamics> nmpc_argmin{dynamics, config};

    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;

    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(horizon) + " ws=" + ws_label;

    // Timing: NLopt
    bench.warmup(5).minEpochIterations(5).title(title)
        .run("nlopt_isres",
             [&]
             {
                 auto u = nmpc_nlopt.solve(x0);
                 ankerl::nanobench::doNotOptimizeAway(u);
             });

    // Timing: argmin
    bench.run("argmin_isres",
              [&]
              {
                  auto u = nmpc_argmin.solve(x0);
                  ankerl::nanobench::doNotOptimizeAway(u);
              });

    // Quality: single solve each for metrics
    ctrlpp::nmpc<double, NX, NU, NloptSolver, Dynamics> q_nlopt{dynamics, config};
    ctrlpp::nmpc<double, NX, NU, ArgminIsres, Dynamics> q_argmin{dynamics, config};

    q_nlopt.solve(x0);
    auto nlopt_diag = q_nlopt.diagnostics();
    auto nlopt_grad = compute_gradient_norm<double, NX, NU>(q_nlopt);

    q_argmin.solve(x0);
    auto argmin_diag = q_argmin.diagnostics();
    auto argmin_grad = compute_gradient_norm<double, NX, NU>(q_argmin);

    quality_metrics nlopt_qm{
        .objective = nlopt_diag.cost,
        .max_constraint_violation = nlopt_diag.max_constraint_violation,
        .gradient_norm = nlopt_grad,
        .success = (nlopt_diag.status == ctrlpp::solve_status::optimal),
        .iterations = nlopt_diag.iterations,
        .solve_time_ms = nlopt_diag.solve_time * 1000.0,
    };

    quality_metrics argmin_qm{
        .objective = argmin_diag.cost,
        .max_constraint_violation = argmin_diag.max_constraint_violation,
        .gradient_norm = argmin_grad,
        .success = (argmin_diag.status == ctrlpp::solve_status::optimal),
        .iterations = argmin_diag.iterations,
        .solve_time_ms = argmin_diag.solve_time * 1000.0,
    };

    write_quality_csv_row(quality_csv, system_name, "nlopt", "isres", "cold",
                          static_cast<int>(NX), horizon, nlopt_qm);
    write_quality_csv_row(quality_csv, system_name, "argmin", "isres", ws_label,
                          static_cast<int>(NX), horizon, argmin_qm);
}

// ---------------------------------------------------------------------------
// Convergence reliability runner
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_convergence(const std::string& system_name,
                     Dynamics dynamics,
                     int horizon,
                     std::ostream& quality_csv)
{
    auto config = make_nmpc_config<NX, NU>(horizon);
    constexpr int num_trials = 100;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-2.0, 2.0);

    int nlopt_successes = 0;
    int argmin_successes = 0;

    for(int trial = 0; trial < num_trials; ++trial)
    {
        Eigen::Matrix<double, NX, 1> x0;
        for(std::size_t i = 0; i < NX; ++i)
            x0(static_cast<Eigen::Index>(i)) = dist(rng);

        ctrlpp::nlopt_settings<double> nlopt_cfg{};
        nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::isres;
        nlopt_cfg.max_eval = 2000;
        NloptSolver nlopt_s{nlopt_cfg};

        ctrlpp::argmin_settings<double> argmin_cfg{};
        argmin_cfg.max_eval = 2000;
        ArgminIsres argmin_s{argmin_cfg};

        ctrlpp::nmpc<double, NX, NU, NloptSolver, Dynamics> nmpc_nlopt{dynamics, config};
        ctrlpp::nmpc<double, NX, NU, ArgminIsres, Dynamics> nmpc_argmin{dynamics, config};

        auto u_nlopt = nmpc_nlopt.solve(x0);
        if(u_nlopt.has_value())
            ++nlopt_successes;

        auto u_argmin = nmpc_argmin.solve(x0);
        if(u_argmin.has_value())
            ++argmin_successes;
    }

    double nlopt_rate = static_cast<double>(nlopt_successes) / num_trials;
    double argmin_rate = static_cast<double>(argmin_successes) / num_trials;

    quality_metrics nlopt_qm{
        .objective = nlopt_rate,
        .max_constraint_violation = 0.0,
        .gradient_norm = 0.0,
        .success = true,
        .iterations = num_trials,
        .solve_time_ms = 0.0,
    };

    quality_metrics argmin_qm{
        .objective = argmin_rate,
        .max_constraint_violation = 0.0,
        .gradient_norm = 0.0,
        .success = true,
        .iterations = num_trials,
        .solve_time_ms = 0.0,
    };

    write_quality_csv_row(quality_csv, system_name, "nlopt", "isres", "convergence",
                          static_cast<int>(NX), horizon, nlopt_qm);
    write_quality_csv_row(quality_csv, system_name, "argmin", "isres", "convergence",
                          static_cast<int>(NX), horizon, argmin_qm);
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_isres_timing.csv");
    std::ofstream quality_csv("bench_isres_quality.csv");
    write_quality_csv_header(quality_csv);

    // Warm-start sweep (BENCH-04): double_integrator_4 and pendulum_2
    for(auto ws : {ctrlpp::warm_start_mode::cold,
                   ctrlpp::warm_start_mode::primal_only,
                   ctrlpp::warm_start_mode::curvature})
    {
        run_benchmark<4, 2>("double_integrator", double_integrator_4, 10, bench, quality_csv, ws);
        run_benchmark<2, 1>("pendulum", pendulum_2, 10, bench, quality_csv, ws);
    }

    // Size sweep (BENCH-05): NX={2,4} only, horizons {10, 20} only (ISRES is slow)
    for(int h : {10, 20})
    {
        run_benchmark<2, 1>("double_integrator", double_integrator_2, h, bench, quality_csv,
                            ctrlpp::warm_start_mode::cold);
        run_benchmark<4, 2>("double_integrator", double_integrator_4, h, bench, quality_csv,
                            ctrlpp::warm_start_mode::cold);
    }

    // Convergence reliability (BENCH-09): 100 trials per solver
    run_convergence<4, 2>("double_integrator", double_integrator_4, 10, quality_csv);
    run_convergence<2, 1>("pendulum", pendulum_2, 10, quality_csv);

    // Render timing CSV
    bench.render(comma_csv_tpl, timing_csv);
}
