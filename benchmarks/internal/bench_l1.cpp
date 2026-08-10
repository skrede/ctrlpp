#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/l1.h"

#include <iostream>
#include <fstream>

int main(int argc, char** argv)
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
    auto controller = ctrlpp::bench::built_or_exit(
        ctrlpp::l1_controller<double, 1, 1>::create(cfg, cutoff_hz, sample_hz), "controller");

    ctrlpp::Vector<double, 1> x{0.5};
    ctrlpp::Vector<double, 1> r{1.0};

    // A rejected cycle performs a fraction of the work, so a run that included
    // one would report a meaningless figure. Establish outside the measured
    // region that the cycle runs, then feed the result to the optimizer barrier
    // inside it so it is neither discarded nor branched on while the clock is
    // running.
    if(const auto stepped = controller.evaluate(x, r); !stepped)
    {
        std::cerr << "bench_l1: l1_controller::evaluate rejected the cycle; the reported figures would be meaningless\n";
        return 1;
    }

    ankerl::nanobench::Bench bench;
    bench.title("L1")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_single_implementation_row(bench, "l1::evaluate", [&] {
        ankerl::nanobench::doNotOptimizeAway(controller.evaluate(x, r).has_value());
    });

    std::ofstream csv("bench_l1.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
