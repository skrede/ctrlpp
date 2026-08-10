#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "ctrlpp/dsp/fir.h"
#include "ctrlpp/dsp/biquad.h"

#include <cmath>
#include <array>
#include <fstream>

namespace
{

constexpr std::size_t fir_taps = 32;
constexpr int settling_samples = 4000;

// The transfer function at zero frequency is the one point on the unit circle a
// recursive filter reaches exactly in finite time, so the identity carries no
// truncation floor: drive the filter with a constant and its output must settle
// on (b0+b1+b2)/(1+a1+a2) formed from its own coefficients.
double dc_gain_residual(ctrlpp::biquad<double>& filter)
{
    const auto c = filter.coefficients();
    const double analytic = (c.b0 + c.b1 + c.b2) / (1.0 + c.a1 + c.a2);
    filter.reset();
    double settled = 0.0;
    for(int i = 0; i < settling_samples; ++i)
        settled = filter.process(1.0);
    return std::abs(settled - analytic) / std::abs(analytic);
}

// A finite impulse response IS the tap vector, so the identity is exact rather
// than asymptotic and no reference implementation enters it.
double impulse_response_residual(ctrlpp::fir<double, fir_taps>& filter,
                                 const std::array<double, fir_taps>& taps)
{
    double worst = 0.0;
    double scale = 0.0;
    for(std::size_t i = 0; i < fir_taps; ++i)
    {
        worst = std::max(worst, std::abs(filter.process(i == 0 ? 1.0 : 0.0) - taps[i]));
        scale = std::max(scale, std::abs(taps[i]));
    }
    return worst / scale;
}

}

int main(int argc, char** argv)
{
    // Biquad low-pass at 10 Hz, 1 kHz sample rate
    auto bq = ctrlpp::biquad<double>::low_pass(10.0, 1000.0).value();
    double sample = 0.5;

    // 32-tap FIR (simple moving average for benchmarking purposes)
    std::array<double, fir_taps> taps{};
    for(auto& t : taps)
        t = 1.0 / double(fir_taps);
    ctrlpp::fir<double, fir_taps> fir_filter(taps);

    const double biquad_residual = dc_gain_residual(bq);
    const double fir_residual = impulse_response_residual(fir_filter, taps);

    ankerl::nanobench::Bench bench;
    bench.title("DSP")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_single_implementation_row(bench, "biquad::process", [&] {
        auto y = bq.process(sample);
        ankerl::nanobench::doNotOptimizeAway(y);
    });

    ctrlpp::bench::run_single_implementation_row(bench, "fir::process", [&] {
        auto y = fir_filter.process(sample);
        ankerl::nanobench::doNotOptimizeAway(y);
    });

    ctrlpp::bench::run_own_criterion_row(
        bench, "relative deviation of this filter's settled constant response from the zero-frequency gain of its own coefficients",
        "biquad::process", biquad_residual, [&] {
            auto y = bq.process(sample);
            ankerl::nanobench::doNotOptimizeAway(y);
        });

    ctrlpp::bench::run_own_criterion_row(
        bench, "max relative deviation of this filter's own impulse response from its tap vector",
        "fir::process", fir_residual, [&] {
            auto y = fir_filter.process(sample);
            ankerl::nanobench::doNotOptimizeAway(y);
        });

    std::ofstream csv("bench_dsp.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
