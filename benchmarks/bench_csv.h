#ifndef HPP_GUARD_BENCHMARKS_BENCH_CSV_H
#define HPP_GUARD_BENCHMARKS_BENCH_CSV_H

// This header includes nanobench, so a translation unit that emits nanobench's
// implementation must define ANKERL_NANOBENCH_IMPLEMENT and include <nanobench.h>
// above this include; defining the macro afterwards leaves the implementation
// unemitted and the link fails.
#include <nanobench.h>

#include <limits>
#include <string>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <utility>

namespace ctrlpp::bench
{

// nanobench resolves {{context(...)}} through std::unordered_map::at, so a row
// whose accuracy entries were never set aborts the render instead of writing an
// empty cell. A blank there would read as a measurement that was taken and came
// out empty, which is the one thing this schema must not be able to say.
constexpr char const* csv_tpl = R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total","accuracy_metric","accuracy_value"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}},"{{context(accuracy_metric)}}","{{context(accuracy_value)}}"
{{/result}})TEMPLATE";

constexpr char const* single_implementation_metric = "no second implementation to compare against";
constexpr char const* single_implementation_value  = "n/a (single implementation)";
constexpr char const* competitor_absent_metric     = "competitor not installed on the measuring station";
constexpr char const* competitor_absent_value      = "n/a (competitor absent)";

inline void report_accuracy(ankerl::nanobench::Bench& bench, char const* metric, double value)
{
    std::ostringstream text;
    text << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) << value;
    bench.context("accuracy_metric", metric).context("accuracy_value", text.str());
}

inline std::string own_criterion_row(char const* arm_name)
{
    return std::string{arm_name} + " own criterion";
}

template <typename OpA, typename OpB>
void run_own_criterion_pair(ankerl::nanobench::Bench& bench, char const* metric, char const* name_a, double value_a,
                            OpA&& op_a, char const* name_b, double value_b, OpB&& op_b)
{
    report_accuracy(bench, metric, value_a);
    bench.run(own_criterion_row(name_a), std::forward<OpA>(op_a));
    report_accuracy(bench, metric, value_b);
    bench.run(own_criterion_row(name_b), std::forward<OpB>(op_b));
}

inline void report_single_implementation(ankerl::nanobench::Bench& bench)
{
    bench.context("accuracy_metric", single_implementation_metric)
        .context("accuracy_value", single_implementation_value);
}

inline void report_competitor_not_installed(ankerl::nanobench::Bench& bench)
{
    bench.context("accuracy_metric", competitor_absent_metric)
        .context("accuracy_value", competitor_absent_value);
}

inline void apply_smoke_switch(ankerl::nanobench::Bench& bench, int32_t argc, char const* const* argv)
{
    for(int32_t i = 1; i < argc; ++i)
        if(std::strcmp(argv[i], "--smoke") == 0)
            bench.warmup(0).epochs(1).epochIterations(1);
}

}

#endif
