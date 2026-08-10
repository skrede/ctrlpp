#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "ctrlpp/control/mrac.h"

#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
    // SISO MRAC: NX=1, NU=1, reference model pole at 0.9
    ctrlpp::mrac_config<double, 1, 1> cfg{};
    cfg.reference_model.A(0, 0) = 0.9;
    cfg.reference_model.B(0, 0) = 0.1;
    cfg.reference_model.C(0, 0) = 1.0;
    cfg.reference_model.D(0, 0) = 0.0;
    cfg.gamma_x(0, 0) = 100.0;
    cfg.gamma_r(0, 0) = 100.0;
    cfg.sign_b(0, 0) = 1.0;

    ctrlpp::mrac_controller<double, 1, 1> controller(cfg);

    ctrlpp::Vector<double, 1> x{0.5};
    ctrlpp::Vector<double, 1> r{1.0};

    // A rejected cycle performs a fraction of the work, so a run that included
    // one would report a meaningless figure. Establish outside the measured
    // region that the cycle runs, then feed the result to the optimizer barrier
    // inside it so it is neither discarded nor branched on while the clock is
    // running.
    if(const auto stepped = controller.evaluate(x, r); !stepped)
    {
        std::cerr << "bench_mrac: mrac_controller::evaluate rejected the cycle; the reported figures would be meaningless\n";
        return 1;
    }

    ankerl::nanobench::Bench bench;
    bench.title("MRAC")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_single_implementation_row(bench, "mrac::evaluate", [&] {
        ankerl::nanobench::doNotOptimizeAway(controller.evaluate(x, r).has_value());
    });

    std::ofstream csv("bench_mrac.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
