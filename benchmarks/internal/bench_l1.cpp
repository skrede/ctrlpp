#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/l1.h"

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // SISO L1: NX=1, NU=1, predictor model with pole at 0.9
    ctrlpp::l1_config<double, 1, 1> cfg{};
    cfg.predictor_model.A(0, 0) = 0.9;
    cfg.predictor_model.B(0, 0) = 0.1;
    cfg.predictor_model.C(0, 0) = 1.0;
    cfg.predictor_model.D(0, 0) = 0.0;
    cfg.gamma(0, 0) = 5000.0;
    cfg.theta_min(0) = -10.0;
    cfg.theta_max(0) = 10.0;

    constexpr double cutoff_hz = 10.0;
    constexpr double sample_hz = 100.0;
    ctrlpp::l1_controller<double, 1, 1> controller(cfg, cutoff_hz, sample_hz);

    ctrlpp::Vector<double, 1> x{0.5};
    ctrlpp::Vector<double, 1> r{1.0};

    ankerl::nanobench::Bench bench;
    bench.title("L1")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true)
        .run("l1::evaluate", [&] {
            auto u = controller.evaluate(x, r);
            ankerl::nanobench::doNotOptimizeAway(u);
        });

    std::ofstream csv("bench_l1.csv");
    bench.render(csv_tpl, csv);
}
