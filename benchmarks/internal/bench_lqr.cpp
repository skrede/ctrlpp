#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/lqr.h"

#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

int main()
{
    // Double integrator: x1(k+1) = x1(k) + dt*x2(k), x2(k+1) = x2(k) + dt*u(k)
    constexpr double dt = 0.01;
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    Eigen::Vector2d B{0.0, dt};
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R = Eigen::Matrix<double, 1, 1>::Identity();

    // Precompute gain for lqr::compute benchmark
    auto K_opt = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    ctrlpp::lqr<double, 2, 1> controller(*K_opt);
    Eigen::Vector2d x{1.0, 0.5};

    ankerl::nanobench::Bench bench;
    bench.title("LQR")
        .performanceCounters(true);

    bench.warmup(50)
        .minEpochIterations(1000)
        .run("lqr_gain", [&] {
            auto K = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
            ankerl::nanobench::doNotOptimizeAway(K);
        });

    bench.warmup(100)
        .minEpochIterations(10000)
        .run("lqr::compute", [&] {
            auto u = controller.compute(x);
            ankerl::nanobench::doNotOptimizeAway(u);
        });

    std::ofstream csv("bench_lqr.csv");
    bench.render(csv_tpl, csv);
}
