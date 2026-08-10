// Competitive benchmark: ctrlpp::online_planner_3rd vs ruckig::Ruckig
// Problem: 1-DOF jerk-constrained time-optimal point-to-point motion, with profile
// synthesis and profile evaluation timed as separate, separately labeled rows.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "trajectory_scan.h"
#include "ctrlpp_trajectory_arm.h"
#include "ruckig_trajectory_arm.h"

#include "bench_csv.h"

#include <ruckig/ruckig.hpp>

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <fstream>

namespace
{

using ctrlpp::bench::comparison_figures;
using ctrlpp::bench::ctrlpp_trajectory_arm;
using ctrlpp::bench::motion_command;
using ctrlpp::bench::ruckig_trajectory_arm;

constexpr char const* library_synthesis  = "ctrlpp::online_planner_3rd::update (synthesis)";
constexpr char const* rival_synthesis    = "ruckig::Ruckig::calculate (synthesis)";
constexpr char const* library_evaluation = "ctrlpp::online_planner_3rd::sample (evaluation)";
constexpr char const* rival_evaluation   = "ruckig::Trajectory::at_time (evaluation)";

constexpr char const* duration_metric =
    "absolute difference of the two arms' total profile durations";
constexpr char const* boundary_metric =
    "absolute position error of this arm's own profile against the commanded target, at the last "
    "representable instant inside its own reported duration";
constexpr char const* margin_metric =
    "worst relative margin of this arm's own profile against the shared velocity, acceleration and "
    "jerk limits, net of the scan's own resolution (negative certifies admissibility)";
constexpr char const* deviation_metric =
    "max relative deviation of the two arms' sampled position, velocity and acceleration over the "
    "common interval";
constexpr char const* symmetry_metric =
    "max relative violation of the point symmetry this arm's own profile has about its own midpoint";

void prime_or_exit(ctrlpp_trajectory_arm& library, ruckig_trajectory_arm& rival)
{
    if(library.synthesize().has_value() && rival.synthesize() == ruckig::Result::Working)
        return;
    std::fprintf(stderr, "benchmark configuration rejected: %s or %s\n", library_synthesis,
                 rival_synthesis);
    std::exit(EXIT_FAILURE);
}

void configure(ankerl::nanobench::Bench& bench, int32_t argc, char const* const* argv)
{
    bench.title("Trajectory: ctrlpp vs ruckig (synthesis and evaluation timed apart)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);
}

void arm_for_timing(ctrlpp_trajectory_arm& library, ruckig_trajectory_arm& rival, double horizon)
{
    library.rewind();
    library.set_horizon(horizon);
    rival.set_horizon(horizon);
}

void emit_synthesis_rows(ankerl::nanobench::Bench& bench, ctrlpp_trajectory_arm& library,
                         ruckig_trajectory_arm& rival, comparison_figures const& figures,
                         motion_command const& cmd)
{
    auto build_library = [&] { ankerl::nanobench::doNotOptimizeAway(library.synthesize()); };
    auto build_rival = [&] { ankerl::nanobench::doNotOptimizeAway(rival.synthesize()); };

    ctrlpp::bench::report_accuracy(bench, duration_metric, figures.duration_gap);
    bench.run(library_synthesis, build_library).run(rival_synthesis, build_rival);
    ctrlpp::bench::run_own_criterion_pair(bench, boundary_metric, library_synthesis,
                                          figures.library.end_position, build_library,
                                          rival_synthesis, figures.rival.end_position, build_rival);
    ctrlpp::bench::run_certificate_pair(
        bench, margin_metric, library_synthesis,
        ctrlpp::bench::worst_net_margin(figures.library, cmd, figures.step), build_library,
        rival_synthesis, ctrlpp::bench::worst_net_margin(figures.rival, cmd, figures.step),
        build_rival);
}

void emit_evaluation_rows(ankerl::nanobench::Bench& bench, ctrlpp_trajectory_arm& library,
                          ruckig_trajectory_arm& rival, comparison_figures const& figures)
{
    auto step_library = [&] { ankerl::nanobench::doNotOptimizeAway(library.advance()); };
    auto step_rival = [&] { ankerl::nanobench::doNotOptimizeAway(rival.advance()); };

    ctrlpp::bench::report_accuracy(bench, deviation_metric, figures.sampled_gap);
    bench.run(library_evaluation, step_library).run(rival_evaluation, step_rival);
    ctrlpp::bench::run_own_criterion_pair(bench, symmetry_metric, library_evaluation,
                                          figures.library.symmetry, step_library, rival_evaluation,
                                          figures.rival.symmetry, step_rival);
}

}

int main(int argc, char** argv)
{
    constexpr motion_command cmd{
        .q0 = 0.0, .q1 = 1.0, .v_max = 2.0, .a_max = 5.0, .j_max = 10.0, .control_period = 0.001};

    ctrlpp_trajectory_arm library{cmd};
    ruckig_trajectory_arm rival{cmd};
    prime_or_exit(library, rival);

    comparison_figures const figures = ctrlpp::bench::measure(library, rival, cmd);
    ctrlpp::bench::report_scan(library_synthesis, figures.library, cmd, figures.step);
    ctrlpp::bench::report_scan(rival_synthesis, figures.rival, cmd, figures.step);

    ankerl::nanobench::Bench bench;
    configure(bench, argc, argv);
    arm_for_timing(library, rival, figures.horizon);

    emit_synthesis_rows(bench, library, rival, figures, cmd);
    emit_evaluation_rows(bench, library, rival, figures);

    std::ofstream csv("bench_trajectory_vs_ruckig.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
