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
#include <string>
#include <vector>
#include <utility>
#include <type_traits>

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

constexpr double pend_dt = 0.05;
auto pendulum_2 = [](const Eigen::Vector2d& x,
                     const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    constexpr double g = 9.81;
    constexpr double l = 1.0;
    double theta = x(0);
    double omega = x(1);
    double alpha = -g / l * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + pend_dt * omega, omega + pend_dt * alpha};
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
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;
using ArgminNwSqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp>;
using ArgminFilterSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_slsqp>;
using ArgminFilterNwSqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_nw_sqp>;

// ---------------------------------------------------------------------------
// Closed-loop simulation
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Solver, typename Dynamics>
auto run_closed_loop(
    Dynamics dynamics,
    const ctrlpp::nmpc_config<double, NX, NU>& config,
    Eigen::Vector<double, static_cast<int>(NX)> x0,
    int sim_steps,
    Solver solver) -> std::pair<double, double>
{
    ctrlpp::nmpc_dynamic<double, NX, NU, Solver, Dynamics> controller{dynamics, config, std::move(solver)};
    double total_cost = 0.0;
    bool all_success = true;

    for(int k = 0; k < sim_steps; ++k)
    {
        auto u = controller.solve(x0);
        if(!u.has_value())
        {
            all_success = false;
            break;
        }
        total_cost += x0.squaredNorm() + u->input.squaredNorm() * 0.1;
        x0 = dynamics(x0, u->input);
    }

    double final_norm = all_success ? x0.norm() : -1.0;
    return {final_norm, total_cost};
}

// ---------------------------------------------------------------------------
// Benchmark runner
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_nmpc_benchmark(const std::string& system_name,
                        Dynamics dynamics,
                        const ctrlpp::nmpc_config<double, NX, NU>& config,
                        const Eigen::Vector<double, static_cast<int>(NX)>& x0,
                        int sim_steps,
                        ankerl::nanobench::Bench& bench,
                        std::ostream& quality_csv)
{
    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(config.horizon)
               + " steps=" + std::to_string(sim_steps);

    auto write_row = [&](char const* label, double norm, double cost)
    {
        quality_csv << system_name << ',' << label << ',' << NX << ',' << config.horizon
                    << ',' << sim_steps << ',' << norm << ',' << cost
                    << ',' << (norm >= 0.0 ? 1 : 0) << '\n';
    };

    bench.title(title);

    // NLopt SLSQP baseline
    {
        ctrlpp::nlopt_settings<double> nlopt_cfg{};
        nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::slsqp;
        bench.run("nlopt_slsqp",
                  [&]
                  {
                      auto result = run_closed_loop<NX, NU>(dynamics, config, x0, sim_steps, NloptSolver{nlopt_cfg});
                      ankerl::nanobench::doNotOptimizeAway(result);
                  });
        auto [norm, cost] = run_closed_loop<NX, NU>(dynamics, config, x0, sim_steps, NloptSolver{nlopt_cfg});
        write_row("nlopt_slsqp", norm, cost);
    }

    // Common argmin settings: bound per-step solve so a non-converging variant
    // cannot stall the closed-loop run. 0.5s per step times sim_steps caps the
    // worst-case wall per single closed-loop iteration.
    ctrlpp::argmin_settings<double> argmin_cfg{};
    argmin_cfg.max_time = 0.5;

    auto bench_argmin_variant = [&]<typename Solver>(std::type_identity<Solver>,
                                                      char const* bench_name,
                                                      bool include_in_bench)
    {
        if (include_in_bench)
        {
            bench.run(bench_name,
                      [&]
                      {
                          auto result = run_closed_loop<NX, NU>(dynamics, config, x0, sim_steps, Solver{argmin_cfg});
                          ankerl::nanobench::doNotOptimizeAway(result);
                      });
        }
        auto [norm, cost] = run_closed_loop<NX, NU>(dynamics, config, x0, sim_steps, Solver{argmin_cfg});
        write_row(bench_name, norm, cost);
    };

    bench_argmin_variant(std::type_identity<ArgminSlsqp>{},        "argmin_slsqp",         true);
    bench_argmin_variant(std::type_identity<ArgminNwSqp>{},        "argmin_nw_sqp",        true);
    bench_argmin_variant(std::type_identity<ArgminFilterSlsqp>{},  "argmin_filter_slsqp",  true);
    // filter_nw_sqp tends to hit max_time on every step on these cells; nanobench
    // batched closed-loop balloons. Single-shot quality only.
    bench_argmin_variant(std::type_identity<ArgminFilterNwSqp>{},  "argmin_filter_nw_sqp", false);
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true).warmup(10).minEpochIterations(5);

    std::ofstream timing_csv("bench_nmpc_timing.csv");
    std::ofstream quality_csv("bench_nmpc_quality.csv");
    quality_csv << "system,solver,nx,horizon,sim_steps,final_state_norm,total_cost,success\n";

    // Double integrator NX=2
    for (int h : {10, 20})
    {
        auto config = make_nmpc_config<2, 1>(h);
        config.Q = Eigen::Matrix2d::Identity() * 10.0;
        config.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
        Eigen::Vector2d x0{1.0, 0.0};
        run_nmpc_benchmark<2, 1>("double_integrator", double_integrator_2, config, x0, 50, bench, quality_csv);
    }

    // Pendulum NX=2
    for (int h : {10, 20})
    {
        auto config = make_nmpc_config<2, 1>(h);
        config.Q = Eigen::Matrix2d::Identity();
        config.R = Eigen::Matrix<double, 1, 1>::Identity();
        Eigen::Vector2d x0{0.3, 0.0};
        run_nmpc_benchmark<2, 1>("pendulum", pendulum_2, config, x0, 100, bench, quality_csv);
    }

    // Double integrator NX=4
    for (int h : {10, 20})
    {
        auto config = make_nmpc_config<4, 2>(h);
        config.Q = Eigen::Matrix4d::Identity();
        config.R = Eigen::Matrix2d::Identity() * 0.1;
        Eigen::Vector4d x0{1.0, 0.0, -0.5, 0.0};
        run_nmpc_benchmark<4, 2>("double_integrator", double_integrator_4, config, x0, 50, bench, quality_csv);
    }

    // Render timing CSV
    bench.render(comma_csv_tpl, timing_csv);
}
