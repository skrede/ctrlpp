#include "bench_metrics.h"
#include "bench_construct.h"

#include "ctrlpp/nmhe.h"
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
// Pendulum dynamics and measurement (NX=2, NU=1, NY=1)
// ---------------------------------------------------------------------------

struct pendulum_dynamics
{
    double dt = 0.05;
    double g = 9.81;
    double l = 1.0;

    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        double theta = x(0);
        double omega = x(1);
        double alpha = -g / l * std::sin(theta) + u(0);
        return Eigen::Vector2d{theta + dt * omega, omega + dt * alpha};
    }
};

struct angle_measurement
{
    auto operator()(const Eigen::Vector2d& x) const -> Eigen::Matrix<double, 1, 1>
    {
        return (Eigen::Matrix<double, 1, 1>() << x(0)).finished();
    }
};

// ---------------------------------------------------------------------------
// Double integrator dynamics and measurement (NX=2, NU=1, NY=1)
// ---------------------------------------------------------------------------

struct double_integrator_dynamics
{
    double dt = 0.1;

    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)};
    }
};

struct position_measurement
{
    auto operator()(const Eigen::Vector2d& x) const -> Eigen::Matrix<double, 1, 1>
    {
        return (Eigen::Matrix<double, 1, 1>() << x(0)).finished();
    }
};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

// ---------------------------------------------------------------------------
// NMHE estimation benchmark
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N,
          typename Solver, typename Dynamics, typename Measurement>
auto run_nmhe_benchmark(
    Dynamics dynamics,
    Measurement measurement,
    const ctrlpp::nmhe_config<double, NX, NU, NY, N>& config,
    int estimation_steps) -> std::pair<double, double>
{
    auto estimator = ctrlpp::bench::built_or_exit(
        ctrlpp::nmhe<double, NX, NU, NY, N, Solver, Dynamics, Measurement>::create(dynamics, measurement, config),
        "nmhe");

    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);

    std::mt19937 rng(42);
    std::normal_distribution<double> process_noise(0.0, 0.001);
    std::normal_distribution<double> meas_noise(0.0, 0.01);

    Eigen::Vector<double, nx> x_true = config.x0;
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    double sum_error = 0.0;
    double max_error = 0.0;
    int count = 0;
    int warmup_steps = estimation_steps / 2;

    for(int k = 0; k < estimation_steps; ++k)
    {
        // Simulate true system with process noise
        x_true = dynamics(x_true, u);
        for(int i = 0; i < nx; ++i)
            x_true(i) += process_noise(rng);

        // Generate noisy measurement
        Eigen::Vector<double, ny> y_true = measurement(x_true);
        Eigen::Vector<double, ny> y_noisy = y_true;
        for(int i = 0; i < ny; ++i)
            y_noisy(i) += meas_noise(rng);

        // Feed to estimator
        estimator.predict(u);
        estimator.update(y_noisy);

        // Compute error after warmup
        if(k >= warmup_steps)
        {
            double err = (x_true - estimator.state()).norm();
            sum_error += err;
            max_error = std::max(max_error, err);
            ++count;
        }
    }

    double mean_error = (count > 0) ? sum_error / count : -1.0;
    return {mean_error, max_error};
}

// ---------------------------------------------------------------------------
// LM policy conditional support
// ---------------------------------------------------------------------------

// LM policy disabled: argmin lm_policy.h has unqualified concept names
// (finite_difference.h bug). Re-enable when argmin fixes this.
#if 0
#include <argmin/solver/lm_policy.h>
#define HAS_LM_POLICY 1

struct argmin_lm_local
{
    using algorithm = argmin::lm_policy<>;
};
using ArgminLm = ctrlpp::argmin_solver<double, argmin_lm_local>;
#endif

// ---------------------------------------------------------------------------
// Benchmark runner helper
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N,
          typename Dynamics, typename Measurement>
void run_estimation_benchmark(const std::string& system_name,
                              Dynamics dynamics,
                              Measurement measurement,
                              const ctrlpp::nmhe_config<double, NX, NU, NY, N>& config,
                              int estimation_steps,
                              ankerl::nanobench::Bench& bench,
                              std::ostream& quality_csv)
{
    auto title = system_name + " NX=" + std::to_string(NX) + " N=" + std::to_string(N)
               + " steps=" + std::to_string(estimation_steps);

    // Timing: NLopt
    bench.title(title)
        .run("nlopt_slsqp",
             [&]
             {
                 auto result = run_nmhe_benchmark<NX, NU, NY, N, NloptSolver>(
                     dynamics, measurement, config, estimation_steps);
                 ankerl::nanobench::doNotOptimizeAway(result);
             });

    // Timing: argmin SLSQP
    bench.run("argmin_slsqp",
              [&]
              {
                  auto result = run_nmhe_benchmark<NX, NU, NY, N, ArgminSlsqp>(
                      dynamics, measurement, config, estimation_steps);
                  ankerl::nanobench::doNotOptimizeAway(result);
              });

    // Quality: single run each
    auto [nlopt_mean, nlopt_max] = run_nmhe_benchmark<NX, NU, NY, N, NloptSolver>(
        dynamics, measurement, config, estimation_steps);
    auto [argmin_mean, argmin_max] = run_nmhe_benchmark<NX, NU, NY, N, ArgminSlsqp>(
        dynamics, measurement, config, estimation_steps);

    quality_csv << system_name << ",nlopt,slsqp," << estimation_steps
                << ',' << nlopt_mean << ',' << nlopt_max
                << ',' << (nlopt_mean >= 0.0 ? 1 : 0) << '\n';
    quality_csv << system_name << ",argmin,slsqp," << estimation_steps
                << ',' << argmin_mean << ',' << argmin_max
                << ',' << (argmin_mean >= 0.0 ? 1 : 0) << '\n';

#ifdef HAS_LM_POLICY
    // Timing: argmin LM
    bench.run("argmin_lm",
              [&]
              {
                  auto result = run_nmhe_benchmark<NX, NU, NY, N, ArgminLm>(
                      dynamics, measurement, config, estimation_steps);
                  ankerl::nanobench::doNotOptimizeAway(result);
              });

    auto [lm_mean, lm_max] = run_nmhe_benchmark<NX, NU, NY, N, ArgminLm>(
        dynamics, measurement, config, estimation_steps);

    quality_csv << system_name << ",argmin,lm," << estimation_steps
                << ',' << lm_mean << ',' << lm_max
                << ',' << (lm_mean >= 0.0 ? 1 : 0) << '\n';
#endif
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true).warmup(10).minEpochIterations(5);

    std::ofstream timing_csv("bench_nmhe_timing.csv");
    std::ofstream quality_csv("bench_nmhe_quality.csv");
    quality_csv << "system,solver,algorithm,estimation_steps,mean_error,max_error,success\n";

    // Pendulum NMHE
    {
        ctrlpp::nmhe_config<double, 2, 1, 1, 5> config;
        config.Q = Eigen::Matrix2d::Identity() * 0.01;
        config.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
        config.P0 = Eigen::Matrix2d::Identity() * 10.0;

        run_estimation_benchmark<2, 1, 1, 5>(
            "pendulum", pendulum_dynamics{}, angle_measurement{},
            config, 100, bench, quality_csv);
    }

    // Double integrator NMHE
    {
        ctrlpp::nmhe_config<double, 2, 1, 1, 5> config;
        config.Q = Eigen::Matrix2d::Identity() * 0.01;
        config.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
        config.P0 = Eigen::Matrix2d::Identity() * 10.0;

        run_estimation_benchmark<2, 1, 1, 5>(
            "double_integrator", double_integrator_dynamics{}, position_measurement{},
            config, 100, bench, quality_csv);
    }

    // Render timing CSV
    bench.render(comma_csv_tpl, timing_csv);
}
