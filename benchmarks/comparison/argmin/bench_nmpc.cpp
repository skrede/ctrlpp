// Closed-loop nonlinear predictive control: the reference nonlinear-programming
// library against argmin's sequential quadratic programming variants, each
// driving the same plant from the same initial state for the same number of
// steps.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "arm_accuracy.h"
#include "bench_metrics.h"
#include "nmpc_arm_probe.h"
#include "bench_construct.h"

#include "nmpc/pendulum.h"
#include "nmpc/double_integrator.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <cstddef>
#include <fstream>
#include <utility>
#include <algorithm>

namespace
{

namespace arms = ctrlpp::bench::argmin_arms;
namespace problems = ctrlpp::bench::problems::nmpc;

using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;
using ArgminNwSqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp>;
using ArgminFilterSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_slsqp>;
using ArgminFilterNwSqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_nw_sqp>;

struct closed_loop_trace
{
    double cost;
    double violation;
    double final_norm;
    Eigen::VectorXd inputs;
};

// Bounding each step's solve keeps a variant that stops converging from
// stalling the whole closed loop.
auto bounded_settings() -> ctrlpp::argmin_settings<double>
{
    ctrlpp::argmin_settings<double> cfg{};
    cfg.max_time = 0.5;
    return cfg;
}

template <typename Action>
void for_each_variant(Action&& action)
{
    ctrlpp::nlopt_settings<double> nlopt_cfg{};
    nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::slsqp;
    const auto cfg = bounded_settings();
    action(NloptSolver{nlopt_cfg}, "nlopt_slsqp", true);
    action(ArgminSlsqp{cfg}, "argmin_slsqp", true);
    action(ArgminNwSqp{cfg}, "argmin_nw_sqp", true);
    action(ArgminFilterSlsqp{cfg}, "argmin_filter_slsqp", true);
    action(ArgminFilterNwSqp{cfg}, "argmin_filter_nw_sqp", false);
}

/// The applied input at every step, the cost the closed loop actually accrued
/// under the configured weights, and the worst the controller's own planned
/// trajectory ever departed from the plant it planned against.
template <std::size_t NX, std::size_t NU, typename Solver, typename Dynamics>
auto run_closed_loop(Dynamics dynamics, const ctrlpp::nmpc_config<double, NX, NU>& config,
                     Eigen::Vector<double, static_cast<int>(NX)> x, int steps, Solver solver)
    -> closed_loop_trace
{
    auto controller = ctrlpp::bench::built_or_exit(
        ctrlpp::nmpc_dynamic<double, NX, NU, Solver, Dynamics>::create(dynamics, config, std::move(solver)),
        "controller");
    Eigen::VectorXd inputs = Eigen::VectorXd::Zero(steps * static_cast<int>(NU));
    double cost = 0.0;
    double violation = 0.0;
    int taken = 0;

    for(int k = 0; k < steps; ++k)
    {
        auto applied = controller.solve(x);
        if(!applied.has_value())
            break;
        violation = std::max(violation, compute_constraint_violation(controller));
        inputs.segment(k * static_cast<int>(NU), static_cast<int>(NU)) = applied->input;
        cost += x.dot(config.Q * x) + applied->input.dot(config.R * applied->input);
        x = dynamics(x, applied->input);
        ++taken;
    }
    return {cost, violation, taken == steps ? x.norm() : -1.0, inputs.head(taken * static_cast<int>(NU))};
}

// nanobench clears its accumulated results whenever the title changes, so the
// cell rides in the row name and the title names the whole benchmark.
auto row_label(char const* algorithm, const std::string& system_name, std::size_t nx, int horizon, int steps)
    -> std::string
{
    return std::string{algorithm} + " " + system_name + " NX=" + std::to_string(nx) + " N="
         + std::to_string(horizon) + " steps=" + std::to_string(steps);
}

void write_closed_loop_row(std::ostream& csv, const std::string& system_name, char const* algorithm,
                           std::size_t nx, int horizon, int steps, const closed_loop_trace& trace)
{
    csv << system_name << ',' << algorithm << ',' << nx << ',' << horizon << ',' << steps << ','
        << trace.final_norm << ',' << trace.cost << ',' << (trace.final_norm >= 0.0 ? 1 : 0) << '\n';
}

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_nmpc_benchmark(const std::string& system_name, Dynamics dynamics,
                        const ctrlpp::nmpc_config<double, NX, NU>& config,
                        const Eigen::Vector<double, static_cast<int>(NX)>& x0, int steps,
                        ankerl::nanobench::Bench& bench, std::ostream& quality_csv)
{
    std::vector<arms::arm_answer> answers;
    for_each_variant(
        [&](auto solver, char const* algorithm, bool benched)
        {
            const auto trace = run_closed_loop<NX, NU>(dynamics, config, x0, steps, std::move(solver));
            write_closed_loop_row(quality_csv, system_name, algorithm, NX, config.horizon, steps, trace);
            if(benched)
                answers.push_back({trace.violation, row_label(algorithm, system_name, NX, config.horizon, steps),
                                   trace.inputs});
        });

    const double spread = arms::solution_spread(answers);
    for_each_variant(
        [&](auto solver, char const* algorithm, bool benched)
        {
            const arms::arm_answer* answer =
                arms::answer_for_label(answers, row_label(algorithm, system_name, NX, config.horizon, steps));
            if(!benched || answer == nullptr)
                return;
            auto drive = [&]
            {
                ankerl::nanobench::doNotOptimizeAway(
                    run_closed_loop<NX, NU>(dynamics, config, x0, steps, solver));
            };
            arms::emit_variant_rows(bench, *answer, spread, drive);
        });
}

void run_cells(ankerl::nanobench::Bench& bench, std::ostream& quality_csv)
{
    for(int horizon : {10, 20})
    {
        auto small = problems::make_nmpc_quadratic_config<2, 1>(horizon);
        small.Q = Eigen::Matrix2d::Identity() * 10.0;
        run_nmpc_benchmark<2, 1>("double_integrator", problems::double_integrator_2, small,
                                 Eigen::Vector2d{1.0, 0.0}, 50, bench, quality_csv);
        run_nmpc_benchmark<4, 2>("double_integrator", problems::double_integrator_4,
                                 problems::make_nmpc_quadratic_config<4, 2>(horizon),
                                 Eigen::Vector4d{1.0, 0.0, -0.5, 0.0}, 50, bench, quality_csv);
        auto swinging = problems::make_nmpc_quadratic_config<2, 1>(horizon);
        swinging.R = Eigen::Matrix<double, 1, 1>::Identity();
        run_nmpc_benchmark<2, 1>("pendulum", problems::pendulum_2, swinging, Eigen::Vector2d{0.3, 0.0}, 100,
                                 bench, quality_csv);
    }
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("NMPC: closed-loop regulation")
        .warmup(10)
        .minEpochIterations(5)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    std::ofstream timing_csv("bench_nmpc_timing.csv");
    std::ofstream quality_csv("bench_nmpc_quality.csv");
    quality_csv << "system,solver,nx,horizon,sim_steps,final_state_norm,total_cost,success\n";

    run_cells(bench, quality_csv);

    bench.render(ctrlpp::bench::csv_tpl, timing_csv);
}
