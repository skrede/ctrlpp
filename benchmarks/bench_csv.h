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

inline std::string certificate_row(char const* arm_name)
{
    return std::string{arm_name} + " certificate";
}

// Every row-running helper below funnels through this one, so the set-then-run
// ordering the schema depends on is written once. A row whose accuracy was set
// before some EARLIER row ran inherits that row's figure silently, with exit 0
// and no diagnostic; only never-set aborts. Nothing downstream can detect the
// inherited case, so no benchmark spells the ordering itself.
template <typename Op>
void run_with_accuracy(ankerl::nanobench::Bench& bench, char const* metric, std::string const& name, double value,
                       Op&& op)
{
    report_accuracy(bench, metric, value);
    bench.run(name, std::forward<Op>(op));
}

template <typename Op>
void run_own_criterion_row(ankerl::nanobench::Bench& bench, char const* metric, char const* arm_name, double value,
                           Op&& op)
{
    run_with_accuracy(bench, metric, own_criterion_row(arm_name), value, std::forward<Op>(op));
}

template <typename Op>
void run_certificate_row(ankerl::nanobench::Bench& bench, char const* metric, char const* arm_name, double value,
                         Op&& op)
{
    run_with_accuracy(bench, metric, certificate_row(arm_name), value, std::forward<Op>(op));
}

template <typename OpA, typename OpB>
void run_own_criterion_pair(ankerl::nanobench::Bench& bench, char const* metric, char const* name_a, double value_a,
                            OpA&& op_a, char const* name_b, double value_b, OpB&& op_b)
{
    run_own_criterion_row(bench, metric, name_a, value_a, std::forward<OpA>(op_a));
    run_own_criterion_row(bench, metric, name_b, value_b, std::forward<OpB>(op_b));
}

template <typename OpA, typename OpB>
void run_certificate_pair(ankerl::nanobench::Bench& bench, char const* metric, char const* name_a, double value_a,
                          OpA&& op_a, char const* name_b, double value_b, OpB&& op_b)
{
    run_certificate_row(bench, metric, name_a, value_a, std::forward<OpA>(op_a));
    run_certificate_row(bench, metric, name_b, value_b, std::forward<OpB>(op_b));
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

template <typename Op>
void run_single_implementation_row(ankerl::nanobench::Bench& bench, std::string const& name, Op&& op)
{
    report_single_implementation(bench);
    bench.run(name, std::forward<Op>(op));
}

template <typename Op>
void run_competitor_absent_row(ankerl::nanobench::Bench& bench, std::string const& name, Op&& op)
{
    report_competitor_not_installed(bench);
    bench.run(name, std::forward<Op>(op));
}

inline void apply_smoke_switch(ankerl::nanobench::Bench& bench, int32_t argc, char const* const* argv)
{
    for(int32_t i = 1; i < argc; ++i)
        if(std::strcmp(argv[i], "--smoke") == 0)
            bench.warmup(0).epochs(1).epochIterations(1);
}

}

#endif
