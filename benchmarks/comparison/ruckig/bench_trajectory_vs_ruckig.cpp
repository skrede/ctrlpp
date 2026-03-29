// Competitive benchmark: ctrlpp::online_planner_3rd vs ruckig::Ruckig
// Problem: 1-DOF jerk-constrained online trajectory generation

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <ruckig/ruckig.hpp>

#include <fstream>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

} // namespace

int main()
{
    // Shared kinematic limits
    constexpr double v_max = 2.0;
    constexpr double a_max = 5.0;
    constexpr double j_max = 10.0;
    constexpr double dt = 0.001;
    constexpr double target = 1.0;

    // ---- ctrlpp setup ----
    ctrlpp::online_planner_3rd<double> ctrlpp_planner({
        .v_max = v_max,
        .a_max = a_max,
        .j_max = j_max,
    });
    ctrlpp_planner.update(target);

    double t_ctrlpp = 0.0;
    // Warm up
    ctrlpp_planner.sample(t_ctrlpp);

    // ---- ruckig setup ----
    ruckig::Ruckig<1> ruckig_otg{dt};
    ruckig::InputParameter<1> input;
    ruckig::OutputParameter<1> output;

    input.current_position = {0.0};
    input.current_velocity = {0.0};
    input.current_acceleration = {0.0};
    input.target_position = {target};
    input.target_velocity = {0.0};
    input.target_acceleration = {0.0};
    input.max_velocity = {v_max};
    input.max_acceleration = {a_max};
    input.max_jerk = {j_max};

    // Warm up
    ruckig_otg.update(input, output);

    // ---- Benchmark ----
    ankerl::nanobench::Bench bench;
    bench.title("Trajectory: ctrlpp vs ruckig")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::online_planner_3rd::sample",
             [&]
             {
                 auto pt = ctrlpp_planner.sample(t_ctrlpp);
                 ankerl::nanobench::doNotOptimizeAway(pt);
                 t_ctrlpp += dt;
             })
        .run("ruckig::Ruckig::update",
             [&]
             {
                 auto r = ruckig_otg.update(input, output);
                 ankerl::nanobench::doNotOptimizeAway(r);
                 output.pass_to_input(input);
             });

    std::ofstream csv("bench_trajectory_vs_ruckig.csv");
    bench.render(comma_csv_tpl, csv);
}
