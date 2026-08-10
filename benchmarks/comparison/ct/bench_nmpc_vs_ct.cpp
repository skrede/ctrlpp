// Competitive benchmark: ctrlpp::nmpc_dynamic vs ct::optcon::NLOptConSolver,
// both driving one damped linear oscillator in receding-horizon closed loop over
// an identical explicit-Euler plant. The runtime-horizon controller is the arm
// here, not the compile-time-horizon one the public nmpc name resolves to.
//
// The two sides run different algorithm families -- a sequential quadratic
// program over the whole horizon against a Gauss-Newton multiple-shooting
// sweep -- so the timing is a claim only where the achieved cost agrees. The
// cost each arm's own loop realizes is therefore published beside its speed,
// under one shared objective neither library computes for itself.
//
// Each timed call builds its controller and runs the whole loop, so repeated
// calls repeat the same cold work instead of warm-starting off the last one.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ct_nloc_arm.h"
#include "nloc_problem.h"

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/nmpc.h"

#include "ctrlpp/mpc/nlopt_solver.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdint>
#include <fstream>

namespace
{

using ctrlpp::bench::ct_nloc_arm;
using ctrlpp::bench::loop_outcome;
using ctrlpp::bench::oscillator_horizon;
using ctrlpp::bench::oscillator_initial_state;
using ctrlpp::bench::oscillator_input;
using ctrlpp::bench::oscillator_input_dim;
using ctrlpp::bench::oscillator_input_weight;
using ctrlpp::bench::oscillator_sample_period;
using ctrlpp::bench::oscillator_sim_steps;
using ctrlpp::bench::oscillator_stage_cost;
using ctrlpp::bench::oscillator_state;
using ctrlpp::bench::oscillator_state_dim;
using ctrlpp::bench::oscillator_state_weight;
using ctrlpp::bench::oscillator_step;

auto oscillator_dynamics = [](const oscillator_state& x, const oscillator_input& u) -> oscillator_state
{
    return oscillator_step(x, u);
};

using ctrlpp_controller = ctrlpp::nmpc_dynamic<double, oscillator_state_dim, oscillator_input_dim,
                                               ctrlpp::nlopt_solver<double>, decltype(oscillator_dynamics)>;

// The stage weights carry the sample period because the shared objective is the
// sampled integral of the running cost, which is what the competitor's own cost
// function already minimizes; the terminal weight is undivided on both sides.
ctrlpp::nmpc_config<double, oscillator_state_dim, oscillator_input_dim> ctrlpp_configuration()
{
    return {.horizon = oscillator_horizon,
            .Q = oscillator_sample_period * oscillator_state_weight::Identity(),
            .R = oscillator_sample_period * oscillator_input_weight::Identity(),
            .Qf = oscillator_state_weight::Identity()};
}

loop_outcome run_ctrlpp_loop()
{
    auto controller = ctrlpp::bench::built_or_exit(
        ctrlpp_controller::create(oscillator_dynamics, ctrlpp_configuration()), "ctrlpp nonlinear predictive control");
    oscillator_state x = oscillator_initial_state();
    loop_outcome out{0.0, 0.0};
    for(int32_t k = 0; k < oscillator_sim_steps; ++k)
    {
        const auto step = ctrlpp::bench::built_or_exit(controller.solve(x), "ctrlpp closed-loop step");
        if(k == 0)
            out.first_input = step.input[0];
        out.cost += oscillator_stage_cost(x, step.input);
        x = oscillator_step(x, step.input);
    }
    return out;
}

loop_outcome run_ct_loop()
{
    ct_nloc_arm arm;
    return arm.run();
}

constexpr char const* deviation_metric = "max abs deviation of the two arms' first applied input";
constexpr char const* cost_metric =
    "cost this arm's own closed loop realizes under the shared sampled objective (lower is better)";

void emit_rows(ankerl::nanobench::Bench& bench, const loop_outcome& ctrlpp_out, const loop_outcome& ct_out)
{
    auto solve_ctrlpp = [&] { ankerl::nanobench::doNotOptimizeAway(run_ctrlpp_loop()); };
    auto solve_ct = [&] { ankerl::nanobench::doNotOptimizeAway(run_ct_loop()); };

    ctrlpp::bench::report_accuracy(bench, deviation_metric,
                                   std::abs(ctrlpp_out.first_input - ct_out.first_input));
    bench.run("ctrlpp::nmpc_dynamic", solve_ctrlpp).run("ct::optcon::NLOptConSolver", solve_ct);
    ctrlpp::bench::run_own_criterion_pair(bench, cost_metric, "ctrlpp::nmpc_dynamic", ctrlpp_out.cost, solve_ctrlpp,
                                          "ct::optcon::NLOptConSolver", ct_out.cost, solve_ct);
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("Nonlinear predictive control: ctrlpp vs ct_optcon (damped oscillator, closed loop)")
        .warmup(2)
        .minEpochIterations(5)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    emit_rows(bench, run_ctrlpp_loop(), run_ct_loop());

    std::ofstream csv("bench_nmpc_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
