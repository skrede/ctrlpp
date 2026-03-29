#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/dare.h"

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // 4-state, 2-input discrete system
    constexpr double dt = 0.01;
    Eigen::Matrix4d A;
    A << 1.0, dt, 0.0, 0.0,
         0.0, 1.0, dt, 0.0,
         0.0, 0.0, 1.0, dt,
         0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<double, 4, 2> B;
    B << 0.0, 0.0,
         dt, 0.0,
         0.0, 0.0,
         0.0, dt;

    Eigen::Matrix4d Q = Eigen::Matrix4d::Identity();
    Eigen::Matrix2d R = Eigen::Matrix2d::Identity();

    ankerl::nanobench::Bench bench;
    bench.title("DARE")
        .warmup(50)
        .minEpochIterations(1000)
        .performanceCounters(true)
        .run("dare_solve", [&] {
            auto P = ctrlpp::dare<double, 4, 2>(A, B, Q, R);
            ankerl::nanobench::doNotOptimizeAway(P);
        });

    std::ofstream csv("bench_dare.csv");
    bench.render(csv_tpl, csv);
}
