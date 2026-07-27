#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_construct.h"

#include "ctrlpp/trajectory/cubic_spline.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <fstream>
#include <vector>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // Cubic spline: 10 knots, natural BC
    ctrlpp::cubic_spline<double>::config spline_cfg{
        .times = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},
        .positions = {0.0, 0.5, 1.0, 0.8, 0.3, 0.0, -0.3, -0.8, -1.0, -0.5},
        .bc = ctrlpp::boundary_condition::natural,
    };
    auto spline = ctrlpp::bench::built_or_exit(
        ctrlpp::cubic_spline<double>::create(spline_cfg), "spline");
    double t_eval = 0.45;

    // Online planner: jerk-limited
    ctrlpp::online_planner_3rd<double>::config planner_cfg{
        .v_max = 2.0,
        .a_max = 5.0,
        .j_max = 10.0,
    };
    auto planner = ctrlpp::bench::built_or_exit(
        ctrlpp::online_planner_3rd<double>::create(planner_cfg), "planner");
    planner.update(1.0);
    double t_sample = 0.01;

    ankerl::nanobench::Bench bench;
    bench.title("Trajectory")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);

    bench.run("cubic_spline::evaluate", [&] {
        auto pt = spline.evaluate(t_eval);
        ankerl::nanobench::doNotOptimizeAway(pt);
    });

    bench.run("online_planner_3rd::sample", [&] {
        auto pt = planner.sample(t_sample);
        ankerl::nanobench::doNotOptimizeAway(pt);
        t_sample += 0.001;
    });

    std::ofstream csv("bench_trajectory.csv");
    bench.render(csv_tpl, csv);
}
