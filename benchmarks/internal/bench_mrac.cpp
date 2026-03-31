#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/mrac.h"

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
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

    ankerl::nanobench::Bench bench;
    bench.title("MRAC")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true)
        .run("mrac::evaluate", [&] {
            auto u = controller.evaluate(x, r);
            ankerl::nanobench::doNotOptimizeAway(u);
        });

    std::ofstream csv("bench_mrac.csv");
    bench.render(csv_tpl, csv);
}
