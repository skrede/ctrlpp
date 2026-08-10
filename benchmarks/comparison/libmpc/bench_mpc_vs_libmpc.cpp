// Competitive benchmark: ctrlpp::mpc vs libmpc++ LMPC
// Problem: linear MPC, double integrator (NX=4, NU=2), horizon N=10, box constraints.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "lmpc/double_integrator.h"

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/state_space.h"

#include <mpc/LMPC.hpp>

#include <Eigen/Dense>

#include <cstdio>
#include <fstream>
#include <algorithm>

namespace
{

namespace problems = ctrlpp::bench::problems::lmpc;

constexpr std::size_t NX = 4;
constexpr std::size_t NU = 2;
constexpr std::size_t NY = NX;
constexpr std::size_t NDU = 0;
constexpr int N = 10;
constexpr double tolerance = 1e-3;
constexpr int iteration_budget = 4000;

constexpr char const* deviation_metric = "max abs deviation of the two arms' planned input sequences";
constexpr char const* feasibility_metric =
    "max abs violation of the plant dynamics and the box bounds by this arm's own planned trajectory";
constexpr char const* ctrlpp_label = "ctrlpp::mpc::solve";
constexpr char const* libmpc_label = "libmpc++::LMPC::step";

using ctrlpp_controller = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;
using libmpc_controller = mpc::LMPC<NX, NU, NDU, NY, N, N>;
using plant = ctrlpp::discrete_state_space<double, NX, NU, NY>;
using settings = ctrlpp::mpc_config<double, NX, NU>;

struct plan
{
    Eigen::MatrixXd states;
    Eigen::MatrixXd inputs;
};

// libmpc++ defaults to a 100-iteration cap with no warm start, so a stock
// pairing races a capped, cold-started arm against one that runs to
// convergence. Both sides get the tolerances, budget, warm start and polish
// ctrlpp::osqp_solver defaults to.
auto matched_parameters() -> mpc::LParameters
{
    mpc::LParameters parameters;
    parameters.eps_abs = tolerance;
    parameters.eps_rel = tolerance;
    parameters.maximum_iteration = iteration_budget;
    parameters.enable_warm_start = true;
    parameters.polish = true;
    return parameters;
}

template <std::size_t ROWS>
auto held_over_horizon(const ctrlpp::Vector<double, ROWS>& bound) -> mpc::mat<ROWS, N>
{
    mpc::mat<ROWS, N> held;
    held.colwise() = bound;
    return held;
}

void configure_libmpc(libmpc_controller& controller, const plant& system, const settings& config)
{
    mpc::mat<NX, NX> transition = system.A;
    mpc::mat<NX, NU> actuation = system.B;
    mpc::mat<NY, NX> observation = mpc::mat<NY, NX>::Identity();
    controller.setStateSpaceModel(transition, actuation, observation);

    mpc::mat<NY, N> output_weight = held_over_horizon<NY>(config.Q.diagonal());
    mpc::mat<NU, N> input_weight = held_over_horizon<NU>(config.R.diagonal());
    mpc::mat<NU, N> rate_weight = mpc::mat<NU, N>::Zero();
    controller.setObjectiveWeights(output_weight, input_weight, rate_weight);

    controller.setStateBounds(held_over_horizon<NX>(*config.x_min), held_over_horizon<NX>(*config.x_max));
    controller.setInputBounds(held_over_horizon<NU>(*config.u_min), held_over_horizon<NU>(*config.u_max));
    controller.setOptimizerParameters(matched_parameters());
}

auto ctrlpp_plan(const ctrlpp_controller& controller) -> plan
{
    const auto trajectory = controller.trajectory();
    plan planned{Eigen::MatrixXd(N + 1, NX), Eigen::MatrixXd(N, NU)};
    for(int k = 0; k <= N; ++k)
        planned.states.row(k) = trajectory->first[static_cast<std::size_t>(k)].transpose();
    for(int k = 0; k < N; ++k)
        planned.inputs.row(k) = trajectory->second[static_cast<std::size_t>(k)].transpose();
    return planned;
}

auto libmpc_plan(libmpc_controller& controller) -> plan
{
    const auto sequence = controller.getOptimalSequence();
    return {sequence.state.topRows(N + 1), sequence.input.topRows(N)};
}

auto box_excess(const Eigen::VectorXd& value, const Eigen::VectorXd& low, const Eigen::VectorXd& high) -> double
{
    return std::max((value - high).cwiseMax(0.0).maxCoeff(), (low - value).cwiseMax(0.0).maxCoeff());
}

// The plan claims to be a trajectory of the plant that respects the box bounds.
// Both halves of that claim are checked on the plan the arm itself returned, so
// neither arm's problem is solved a second time here.
auto plan_violation(const plan& planned, const plant& system, const settings& config) -> double
{
    double worst = 0.0;
    for(int k = 0; k < N; ++k)
    {
        const Eigen::VectorXd defect = planned.states.row(k + 1).transpose()
                                     - system.A * planned.states.row(k).transpose()
                                     - system.B * planned.inputs.row(k).transpose();
        worst = std::max({worst, defect.cwiseAbs().maxCoeff(),
                          box_excess(planned.inputs.row(k).transpose(), *config.u_min, *config.u_max)});
    }
    for(int k = 0; k <= N; ++k)
        worst = std::max(worst, box_excess(planned.states.row(k).transpose(), *config.x_min, *config.x_max));
    return worst;
}

void emit_rows(ankerl::nanobench::Bench& bench, ctrlpp_controller& mine, libmpc_controller& theirs,
               const plant& system, const settings& config)
{
    const Eigen::Vector4d x0 = problems::double_integrator_4_2_x0_default();
    const mpc::cvec<NX> their_x0 = x0;
    auto solve_mine = [&] { ankerl::nanobench::doNotOptimizeAway(mine.solve(x0)); };
    auto solve_theirs = [&] { ankerl::nanobench::doNotOptimizeAway(theirs.optimize(their_x0, mpc::cvec<NU>::Zero())); };

    const plan my_plan = ctrlpp_plan(mine);
    const plan their_plan = libmpc_plan(theirs);
    ctrlpp::bench::report_accuracy(bench, deviation_metric,
                                   (my_plan.inputs - their_plan.inputs).cwiseAbs().maxCoeff());
    bench.run(ctrlpp_label, solve_mine).run(libmpc_label, solve_theirs);
    ctrlpp::bench::run_own_criterion_pair(bench, feasibility_metric, ctrlpp_label,
                                          plan_violation(my_plan, system, config), solve_mine, libmpc_label,
                                          plan_violation(their_plan, system, config), solve_theirs);
}

// Without an explicit terminal weight ctrlpp substitutes the infinite-horizon
// Riccati solution, which libmpc++ has no way to express: its output weight is
// one diagonal per step. The two arms would then minimize different costs.
auto matched_config() -> settings
{
    settings config = problems::make_double_integrator_4_2_config(N);
    config.Qf = config.Q;
    return config;
}

}

int main(int argc, char** argv)
{
    const plant system = problems::make_double_integrator_4_2_state_space();
    const settings config = matched_config();

    auto mine = ctrlpp::bench::built_or_exit(ctrlpp_controller::create(system, config), "ctrlpp_mpc");
    if(!mine.solve(problems::double_integrator_4_2_x0_default()).has_value())
    {
        std::fprintf(stderr, "ctrlpp::mpc declined the problem; the reported figures would be meaningless\n");
        return 1;
    }

    libmpc_controller theirs;
    configure_libmpc(theirs, system, config);
    theirs.optimize(problems::double_integrator_4_2_x0_default(), mpc::cvec<NU>::Zero());

    ankerl::nanobench::Bench bench;
    bench.title("MPC: ctrlpp vs libmpc++").warmup(50).minEpochIterations(100).performanceCounters(true).relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);
    emit_rows(bench, mine, theirs, system, config);

    std::ofstream csv("bench_mpc_vs_libmpc.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
