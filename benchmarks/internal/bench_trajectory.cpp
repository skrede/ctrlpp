#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "comparison/ruckig/ctrlpp_trajectory_arm.h"

#include "ctrlpp/trajectory/cubic_spline.h"

#include <cmath>
#include <fstream>
#include <vector>

namespace
{

// The interpolation conditions are what the spline is defined by, so the
// distance between where it evaluates at a knot time and the position that knot
// carries is its own arithmetic and nothing else.
double knot_interpolation_residual(const ctrlpp::cubic_spline<double>& spline,
                                   const ctrlpp::cubic_spline<double>::config& cfg)
{
    double worst = 0.0;
    double span = 0.0;
    for(std::size_t i = 0; i < cfg.times.size(); ++i)
    {
        worst = std::max(worst, std::abs(spline.evaluate(cfg.times[i]).position[0] - cfg.positions[i]));
        span = std::max(span, std::abs(cfg.positions[i]));
    }
    return worst / span;
}

}

int main(int argc, char** argv)
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
    const ctrlpp::bench::motion_command command{
        .q0 = 0.0, .q1 = 1.0, .v_max = 2.0, .a_max = 5.0, .j_max = 10.0, .control_period = 0.001};
    ctrlpp::bench::ctrlpp_trajectory_arm planner(command);
    if(!planner.synthesize().has_value())
        return 1;

    const double spline_residual = knot_interpolation_residual(spline, spline_cfg);
    const double planner_residual =
        ctrlpp::bench::scan_arm(planner, command, ctrlpp::bench::scan_step(command)).symmetry;
    planner.set_horizon(planner.duration());

    ankerl::nanobench::Bench bench;
    bench.title("Trajectory")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_single_implementation_row(bench, "cubic_spline::evaluate", [&] {
        auto pt = spline.evaluate(t_eval);
        ankerl::nanobench::doNotOptimizeAway(pt);
    });

    ctrlpp::bench::run_single_implementation_row(bench, "online_planner_3rd::sample", [&] {
        auto pt = planner.advance();
        ankerl::nanobench::doNotOptimizeAway(pt);
    });

    ctrlpp::bench::run_own_criterion_row(
        bench, "max relative deviation of this spline from the knot positions it interpolates",
        "cubic_spline::evaluate", spline_residual, [&] {
            auto pt = spline.evaluate(t_eval);
            ankerl::nanobench::doNotOptimizeAway(pt);
        });

    ctrlpp::bench::run_own_criterion_row(
        bench, "max relative violation of the point symmetry this arm's own profile has about its own midpoint",
        "online_planner_3rd::sample", planner_residual, [&] {
            auto pt = planner.advance();
            ankerl::nanobench::doNotOptimizeAway(pt);
        });

    std::ofstream csv("bench_trajectory.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
