#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/sysid/n4sid.h"
#include "ctrlpp/sysid/batch_arx.h"

#include <Eigen/Dense>

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // Generate synthetic data from a known 2nd-order discrete system:
    //   y(k) = 1.5*y(k-1) - 0.7*y(k-2) + 0.5*u(k-1) + 0.2*u(k-2)
    constexpr int N = 200;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);
    Y.setZero();
    U.setZero();

    // Step input after sample 5
    for (int k = 0; k < N; ++k)
        U(k) = (k >= 5) ? 1.0 : 0.0;

    for (int k = 2; k < N; ++k)
        Y(k) = 1.5 * Y(k - 1) - 0.7 * Y(k - 2) + 0.5 * U(k - 1) + 0.2 * U(k - 2);

    ankerl::nanobench::Bench bench;
    bench.title("SysID")
        .warmup(10)
        .minEpochIterations(100)
        .performanceCounters(true);

    bench.run("batch_arx::identify", [&] {
        auto result = ctrlpp::batch_arx<2, 2>(Y, U);
        ankerl::nanobench::doNotOptimizeAway(result);
    });

    bench.run("n4sid::identify", [&] {
        auto result = ctrlpp::n4sid<2>(Y, U);
        ankerl::nanobench::doNotOptimizeAway(result);
    });

    std::ofstream csv("bench_sysid.csv");
    bench.render(csv_tpl, csv);
}
