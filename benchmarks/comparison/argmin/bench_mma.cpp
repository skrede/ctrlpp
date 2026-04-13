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
#include <iostream>
#include <random>
#include <stdexcept>
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

constexpr double di8_dt = 0.1;
using Vec8 = Eigen::Matrix<double, 8, 1>;
using Vec4 = Eigen::Vector4d;

auto double_integrator_8 = [](const Vec8& x, const Vec4& u) -> Vec8
{
    Vec8 xn;
    xn(0) = x(0) + di8_dt * x(1);
    xn(1) = x(1) + di8_dt * u(0);
    xn(2) = x(2) + di8_dt * x(3);
    xn(3) = x(3) + di8_dt * u(1);
    xn(4) = x(4) + di8_dt * x(5);
    xn(5) = x(5) + di8_dt * u(2);
    xn(6) = x(6) + di8_dt * x(7);
    xn(7) = x(7) + di8_dt * u(3);
    return xn;
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
using ArgminMma = ctrlpp::argmin_solver<double, ctrlpp::argmin_mma>;
using ArgminGcmma = ctrlpp::argmin_solver<double, ctrlpp::argmin_gcmma>;

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
// Write incompatible row to quality CSV
// ---------------------------------------------------------------------------

void write_incompatible_row(std::ostream& quality_csv,
                            std::string_view system,
                            std::string_view solver,
                            std::string_view warm_start,
                            int nx,
                            int horizon)
{
    quality_metrics qm{
        .objective = 0.0,
        .max_constraint_violation = 0.0,
        .gradient_norm = 0.0,
        .success = false,
        .iterations = 0,
        .solve_time_ms = 0.0,
    };
    write_quality_csv_row(quality_csv, system, solver, "mma_incompatible", warm_start, nx, horizon, qm);
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

    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;

    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(horizon) + " ws=" + ws_label;

    // NLopt MMA: attempt NMPC construction -- will throw on equality constraints
    try
    {
        ctrlpp::nlopt_settings<double> nlopt_cfg{};
        nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::mma;
        NloptSolver nlopt_solver{nlopt_cfg};

        ctrlpp::nmpc<double, NX, NU, NloptSolver, Dynamics> nmpc_nlopt{dynamics, config};
        nmpc_nlopt.solve(x0);

        // If we get here, no equality constraints -- should not happen with NMPC
        auto nlopt_diag = nmpc_nlopt.diagnostics();
        auto nlopt_grad = compute_gradient_norm<double, NX, NU>(nmpc_nlopt);
        quality_metrics nlopt_qm{
            .objective = nlopt_diag.cost,
            .max_constraint_violation = nlopt_diag.max_constraint_violation,
            .gradient_norm = nlopt_grad,
            .success = (nlopt_diag.status == ctrlpp::solve_status::optimal),
            .iterations = nlopt_diag.iterations,
            .solve_time_ms = nlopt_diag.solve_time * 1000.0,
        };
        write_quality_csv_row(quality_csv, system_name, "nlopt", "mma", "cold",
                              static_cast<int>(NX), horizon, nlopt_qm);
    }
    catch(const std::invalid_argument&)
    {
        // Expected: MMA cannot handle NMPC equality constraints
        write_incompatible_row(quality_csv, system_name, "nlopt", "cold",
                               static_cast<int>(NX), horizon);
    }

    // Argmin MMA: attempt NMPC construction -- will throw on equality constraints
    try
    {
        ctrlpp::argmin_settings<double> argmin_cfg{};
        argmin_cfg.warm_start = ws_mode;
        ArgminMma argmin_solver{argmin_cfg};

        ctrlpp::nmpc<double, NX, NU, ArgminMma, Dynamics> nmpc_argmin{dynamics, config};
        nmpc_argmin.solve(x0);

        auto argmin_diag = nmpc_argmin.diagnostics();
        auto argmin_grad = compute_gradient_norm<double, NX, NU>(nmpc_argmin);
        write_quality_csv_row(quality_csv, system_name, "argmin", "mma", ws_label,
                              static_cast<int>(NX), horizon, quality_metrics{
                                  .objective = argmin_diag.cost,
                                  .max_constraint_violation = argmin_diag.max_constraint_violation,
                                  .gradient_norm = argmin_grad,
                                  .success = (argmin_diag.status == ctrlpp::solve_status::optimal),
                                  .iterations = argmin_diag.iterations,
                                  .solve_time_ms = argmin_diag.solve_time * 1000.0,
                              });
    }
    catch(const std::invalid_argument&)
    {
        write_incompatible_row(quality_csv, system_name, "argmin_mma", ws_label,
                               static_cast<int>(NX), horizon);
    }

    // Argmin GCMMA: same constraint limitation as MMA
    try
    {
        ctrlpp::argmin_settings<double> gcmma_cfg{};
        gcmma_cfg.warm_start = ws_mode;
        ArgminGcmma gcmma_solver{gcmma_cfg};

        ctrlpp::nmpc<double, NX, NU, ArgminGcmma, Dynamics> nmpc_gcmma{dynamics, config};
        nmpc_gcmma.solve(x0);

        auto gcmma_diag = nmpc_gcmma.diagnostics();
        auto gcmma_grad = compute_gradient_norm<double, NX, NU>(nmpc_gcmma);
        write_quality_csv_row(quality_csv, system_name, "argmin", "gcmma", ws_label,
                              static_cast<int>(NX), horizon, quality_metrics{
                                  .objective = gcmma_diag.cost,
                                  .max_constraint_violation = gcmma_diag.max_constraint_violation,
                                  .gradient_norm = gcmma_grad,
                                  .success = (gcmma_diag.status == ctrlpp::solve_status::optimal),
                                  .iterations = gcmma_diag.iterations,
                                  .solve_time_ms = gcmma_diag.solve_time * 1000.0,
                              });
    }
    catch(const std::invalid_argument&)
    {
        write_incompatible_row(quality_csv, system_name, "argmin_gcmma", ws_label,
                               static_cast<int>(NX), horizon);
    }
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

    // NLopt MMA: always incompatible with NMPC
    write_incompatible_row(quality_csv, system_name, "nlopt", "convergence",
                           static_cast<int>(NX), horizon);

    // Argmin MMA and GCMMA: incompatible with NMPC equality constraints
    write_incompatible_row(quality_csv, system_name, "argmin_mma", "convergence",
                           static_cast<int>(NX), horizon);
    write_incompatible_row(quality_csv, system_name, "argmin_gcmma", "convergence",
                           static_cast<int>(NX), horizon);
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_mma_timing.csv");
    std::ofstream quality_csv("bench_mma_quality.csv");
    write_quality_csv_header(quality_csv);

    // Warm-start sweep (BENCH-04): only argmin MMA (NLopt MMA is incompatible)
    for(auto ws : {ctrlpp::warm_start_mode::cold,
                   ctrlpp::warm_start_mode::primal_only,
                   ctrlpp::warm_start_mode::curvature})
    {
        run_benchmark<4, 2>("double_integrator", double_integrator_4, 10, bench, quality_csv, ws);
        run_benchmark<2, 1>("pendulum", pendulum_2, 10, bench, quality_csv, ws);
    }

    // Size sweep (BENCH-05): horizons 10, 20, 30 with cold start
    for(int h : {10, 20, 30})
    {
        run_benchmark<2, 1>("double_integrator", double_integrator_2, h, bench, quality_csv,
                            ctrlpp::warm_start_mode::cold);
        run_benchmark<4, 2>("double_integrator", double_integrator_4, h, bench, quality_csv,
                            ctrlpp::warm_start_mode::cold);
        run_benchmark<8, 4>("double_integrator", double_integrator_8, h, bench, quality_csv,
                            ctrlpp::warm_start_mode::cold);
    }

    // Convergence reliability (BENCH-09)
    run_convergence<4, 2>("double_integrator", double_integrator_4, 10, quality_csv);
    run_convergence<2, 1>("pendulum", pendulum_2, 10, quality_csv);

    // Render timing CSV
    bench.render(comma_csv_tpl, timing_csv);
}
