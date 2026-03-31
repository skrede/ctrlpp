#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/pid.h"

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid ctrl(cfg);
    auto sp = Vec::Constant(1.0);
    auto meas = Vec::Constant(0.5);
    constexpr double dt = 0.01;

    ankerl::nanobench::Bench bench;
    bench.title("PID")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true)
        .run("pid::compute", [&] {
            auto u = ctrl.compute(sp, meas, dt);
            ankerl::nanobench::doNotOptimizeAway(u);
        });

    std::ofstream csv("bench_pid.csv");
    bench.render(csv_tpl, csv);
}
