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

// ---------------------------------------------------------------------------
// Closed-loop simulation
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Solver, typename Dynamics>
auto run_closed_loop(
    Dynamics dynamics,
    const ctrlpp::nmpc_config<double, NX, NU>& config,
    Eigen::Vector<double, static_cast<int>(NX)> x0,
    int sim_steps) -> std::pair<double, double>
{
    ctrlpp::nmpc<double, NX, NU, Solver, Dynamics> controller{dynamics, config};
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
        total_cost += x0.squaredNorm() + u->squaredNorm() * 0.1;
        x0 = dynamics(x0, *u);
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

    // Timing: NLopt closed-loop
    bench.title(title)
        .run("nlopt_slsqp",
             [&]
             {
                 auto result = run_closed_loop<NX, NU, NloptSolver>(dynamics, config, x0, sim_steps);
                 ankerl::nanobench::doNotOptimizeAway(result);
             });

    // Timing: argmin closed-loop
    bench.run("argmin_slsqp",
              [&]
              {
                  auto result = run_closed_loop<NX, NU, ArgminSlsqp>(dynamics, config, x0, sim_steps);
                  ankerl::nanobench::doNotOptimizeAway(result);
              });

    // Quality: single run each
    auto [nlopt_norm, nlopt_cost] = run_closed_loop<NX, NU, NloptSolver>(dynamics, config, x0, sim_steps);
    auto [argmin_norm, argmin_cost] = run_closed_loop<NX, NU, ArgminSlsqp>(dynamics, config, x0, sim_steps);

    quality_csv << system_name << ",nlopt_slsqp," << NX << ',' << config.horizon
                << ',' << sim_steps << ',' << nlopt_norm << ',' << nlopt_cost
                << ',' << (nlopt_norm >= 0.0 ? 1 : 0) << '\n';
    quality_csv << system_name << ",argmin_slsqp," << NX << ',' << config.horizon
                << ',' << sim_steps << ',' << argmin_norm << ',' << argmin_cost
                << ',' << (argmin_norm >= 0.0 ? 1 : 0) << '\n';
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
    {
        auto config = make_nmpc_config<2, 1>(10);
        config.Q = Eigen::Matrix2d::Identity() * 10.0;
        config.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
        Eigen::Vector2d x0{1.0, 0.0};
        run_nmpc_benchmark<2, 1>("double_integrator", double_integrator_2, config, x0, 50, bench, quality_csv);
    }

    // Pendulum NX=2
    {
        auto config = make_nmpc_config<2, 1>(10);
        config.Q = Eigen::Matrix2d::Identity();
        config.R = Eigen::Matrix<double, 1, 1>::Identity();
        Eigen::Vector2d x0{0.3, 0.0};
        run_nmpc_benchmark<2, 1>("pendulum", pendulum_2, config, x0, 100, bench, quality_csv);
    }

    // Double integrator NX=4
    {
        auto config = make_nmpc_config<4, 2>(10);
        config.Q = Eigen::Matrix4d::Identity();
        config.R = Eigen::Matrix2d::Identity() * 0.1;
        Eigen::Vector4d x0{1.0, 0.0, -0.5, 0.0};
        run_nmpc_benchmark<4, 2>("double_integrator", double_integrator_4, config, x0, 50, bench, quality_csv);
    }

    // Render timing CSV
    bench.render(comma_csv_tpl, timing_csv);
}
