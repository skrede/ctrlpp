#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/lie/so3.h"

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // Rotation vector: 45 degrees about [1,1,1]/sqrt(3)
    ctrlpp::Vector<double, 3> omega;
    omega << 0.4534498, 0.4534498, 0.4534498;

    // Unit quaternion from that rotation vector (for log benchmark)
    auto q = ctrlpp::so3::exp(omega);

    ankerl::nanobench::Bench bench;
    bench.title("SO3")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);

    bench.run("so3::exp", [&] {
        auto result = ctrlpp::so3::exp(omega);
        ankerl::nanobench::doNotOptimizeAway(result);
    });

    bench.run("so3::log", [&] {
        auto result = ctrlpp::so3::log(q);
        ankerl::nanobench::doNotOptimizeAway(result);
    });

    std::ofstream csv("bench_so3.csv");
    bench.render(csv_tpl, csv);
}
