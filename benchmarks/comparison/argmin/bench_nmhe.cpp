// Nonlinear moving-horizon estimation: the reference nonlinear-programming
// library against argmin's sequential quadratic program, both fed the identical
// deterministic noisy record.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "arm_accuracy.h"
#include "bench_construct.h"

#include "nmpc/pendulum.h"
#include "nmpc/double_integrator.h"

#include "ctrlpp/nmhe.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <random>
#include <cstddef>
#include <fstream>
#include <algorithm>

namespace
{

namespace arms = ctrlpp::bench::argmin_arms;
namespace problems = ctrlpp::bench::problems::nmpc;

using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' estimate sequences";
constexpr char const* error_metric =
    "max deviation of this arm's own estimate from the simulated true state, after warm-up";
constexpr char const* nlopt_label = "nlopt_slsqp";
constexpr char const* argmin_label = "argmin_slsqp";

struct pendulum_dynamics
{
    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        return problems::pendulum_2(x, u);
    }
};

struct double_integrator_dynamics
{
    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        return problems::double_integrator_2(x, u);
    }
};

struct first_state_measurement
{
    auto operator()(const Eigen::Vector2d& x) const -> Eigen::Matrix<double, 1, 1>
    {
        return (Eigen::Matrix<double, 1, 1>() << x(0)).finished();
    }
};

struct estimation_trace
{
    double max_error;
    double mean_error;
    Eigen::VectorXd estimates;
};

/// The record is simulated, so the true state is known exactly at every step and
/// the distance from it is a real error magnitude rather than a fit statistic.
/// The first half of the run is the estimator's transient and is excluded.
template <std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, typename Solver, typename Dynamics,
          typename Measurement>
auto run_estimation(Dynamics dynamics, Measurement measurement,
                    const ctrlpp::nmhe_config<double, NX, NU, NY, N>& config, int steps) -> estimation_trace
{
    auto estimator = ctrlpp::bench::built_or_exit(
        ctrlpp::nmhe<double, NX, NU, NY, N, Solver, Dynamics, Measurement>::create(dynamics, measurement, config),
        "nmhe");
    std::mt19937 rng(42);
    std::normal_distribution<double> process_noise(0.0, 0.001);
    std::normal_distribution<double> measurement_noise(0.0, 0.01);
    Eigen::Vector<double, static_cast<int>(NX)> truth = config.x0;
    const Eigen::Matrix<double, 1, 1> input = Eigen::Matrix<double, 1, 1>::Zero();
    std::vector<double> estimates;
    double sum_error = 0.0;
    double max_error = 0.0;
    int scored = 0;

    for(int k = 0; k < steps; ++k)
    {
        truth = dynamics(truth, input);
        for(int i = 0; i < static_cast<int>(NX); ++i)
            truth(i) += process_noise(rng);
        Eigen::Vector<double, static_cast<int>(NY)> reading = measurement(truth);
        for(int i = 0; i < static_cast<int>(NY); ++i)
            reading(i) += measurement_noise(rng);

        estimator.predict(input);
        // A refused measurement produces no estimate, so a step recorded for it
        // would be one that never happened.
        if(!estimator.update(reading))
            continue;
        if(k < steps / 2)
            continue;

        const double error = (truth - estimator.state()).norm();
        sum_error += error;
        max_error = std::max(max_error, error);
        ++scored;
        for(int i = 0; i < static_cast<int>(NX); ++i)
            estimates.push_back(estimator.state()(i));
    }

    return {max_error, scored > 0 ? sum_error / scored : -1.0,
            Eigen::Map<Eigen::VectorXd>(estimates.data(), static_cast<Eigen::Index>(estimates.size()))};
}

// nanobench clears its accumulated results whenever the title changes, so the
// system rides in the row name and the title names the whole benchmark.
auto row_label(char const* algorithm, const std::string& system_name) -> std::string
{
    return std::string{algorithm} + " " + system_name;
}

void write_estimation_row(std::ostream& csv, const std::string& system_name, char const* algorithm, int steps,
                          const estimation_trace& trace)
{
    csv << system_name << ',' << algorithm << ",slsqp," << steps << ',' << trace.mean_error << ','
        << trace.max_error << ',' << (trace.mean_error >= 0.0 ? 1 : 0) << '\n';
}

template <std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, typename Dynamics, typename Measurement>
void run_estimation_benchmark(const std::string& system_name, Dynamics dynamics, Measurement measurement,
                              const ctrlpp::nmhe_config<double, NX, NU, NY, N>& config, int steps,
                              ankerl::nanobench::Bench& bench, std::ostream& quality_csv)
{
    const auto reference = run_estimation<NX, NU, NY, N, NloptSolver>(dynamics, measurement, config, steps);
    const auto candidate = run_estimation<NX, NU, NY, N, ArgminSlsqp>(dynamics, measurement, config, steps);
    write_estimation_row(quality_csv, system_name, nlopt_label, steps, reference);
    write_estimation_row(quality_csv, system_name, argmin_label, steps, candidate);

    const std::string nlopt_row = row_label(nlopt_label, system_name);
    const std::string argmin_row = row_label(argmin_label, system_name);
    const double deviation = arms::solution_spread(
        {{reference.max_error, nlopt_row, reference.estimates}, {candidate.max_error, argmin_row, candidate.estimates}});

    auto run_reference = [&] {
        ankerl::nanobench::doNotOptimizeAway(
            run_estimation<NX, NU, NY, N, NloptSolver>(dynamics, measurement, config, steps));
    };
    auto run_candidate = [&] {
        ankerl::nanobench::doNotOptimizeAway(
            run_estimation<NX, NU, NY, N, ArgminSlsqp>(dynamics, measurement, config, steps));
    };

    ctrlpp::bench::run_with_accuracy(bench, deviation_metric, nlopt_row, deviation, run_reference);
    ctrlpp::bench::run_with_accuracy(bench, deviation_metric, argmin_row, deviation, run_candidate);
    ctrlpp::bench::run_own_criterion_pair(bench, error_metric, nlopt_row.c_str(), reference.max_error,
                                          run_reference, argmin_row.c_str(), candidate.max_error, run_candidate);
}

auto make_estimation_config() -> ctrlpp::nmhe_config<double, 2, 1, 1, 5>
{
    ctrlpp::nmhe_config<double, 2, 1, 1, 5> config;
    config.Q = Eigen::Matrix2d::Identity() * 0.01;
    config.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
    config.P0 = Eigen::Matrix2d::Identity() * 10.0;
    return config;
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("NMHE: reference SLSQP vs argmin SLSQP")
        .warmup(10)
        .minEpochIterations(5)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    std::ofstream timing_csv("bench_nmhe_timing.csv");
    std::ofstream quality_csv("bench_nmhe_quality.csv");
    quality_csv << "system,solver,algorithm,estimation_steps,mean_error,max_error,success\n";

    const auto config = make_estimation_config();
    run_estimation_benchmark<2, 1, 1, 5>("pendulum", pendulum_dynamics{}, first_state_measurement{}, config, 100,
                                         bench, quality_csv);
    run_estimation_benchmark<2, 1, 1, 5>("double_integrator", double_integrator_dynamics{},
                                         first_state_measurement{}, config, 100, bench, quality_csv);

    bench.render(ctrlpp::bench::csv_tpl, timing_csv);
}
