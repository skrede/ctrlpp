#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/dsp/fir.h"
#include "ctrlpp/dsp/biquad.h"

#include <array>
#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // Biquad low-pass at 10 Hz, 1 kHz sample rate
    auto bq = ctrlpp::biquad<double>::low_pass(10.0, 1000.0).value();
    double sample = 0.5;

    // 32-tap FIR (simple moving average for benchmarking purposes)
    std::array<double, 32> taps{};
    for (auto& t : taps)
        t = 1.0 / 32.0;
    ctrlpp::fir<double, 32> fir_filter(taps);

    ankerl::nanobench::Bench bench;
    bench.title("DSP")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);

    bench.run("biquad::process", [&] {
        auto y = bq.process(sample);
        ankerl::nanobench::doNotOptimizeAway(y);
    });

    bench.run("fir::process", [&] {
        auto y = fir_filter.process(sample);
        ankerl::nanobench::doNotOptimizeAway(y);
    });

    std::ofstream csv("bench_dsp.csv");
    bench.render(csv_tpl, csv);
}
