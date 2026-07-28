#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/pid.h"

#include <iostream>
#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    using Pid = ctrlpp::pid<double, 1>;
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

    // A rejected cycle performs a fraction of the work, so a run that included
    // one would report a meaningless figure. Establish outside the measured
    // region that the cycle runs, then feed the result to the optimizer barrier
    // inside it so it is neither discarded nor branched on while the clock is
    // running.
    if(const auto stepped = ctrl.compute(sp, meas, dt); !stepped)
    {
        std::cerr << "bench_pid: pid::compute rejected the cycle; the reported figures would be meaningless\n";
        return 1;
    }

    ankerl::nanobench::Bench bench;
    bench.title("PID")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true)
        .run("pid::compute", [&] {
            ankerl::nanobench::doNotOptimizeAway(ctrl.compute(sp, meas, dt).has_value());
        });

    std::ofstream csv("bench_pid.csv");
    bench.render(csv_tpl, csv);
}
