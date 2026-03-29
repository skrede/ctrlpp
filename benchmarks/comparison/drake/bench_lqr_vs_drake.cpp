// Competitive benchmark: ctrlpp::lqr_gain vs Drake LinearQuadraticRegulator
// Problem: LQR gain computation for double-integrator (NX=4, NU=2)
//
// Drake requires manual installation. No AUR package. Bazel-only build system.
// Pre-built tar.gz available for Ubuntu/macOS. See benchmarks/README.md.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/lqr.h"

#include <drake/systems/controllers/linear_quadratic_regulator.h>

#include <Eigen/Dense>

#include <fstream>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

} // namespace

int main()
{
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;
    constexpr double dt = 0.1;

    Eigen::Matrix4d A = Eigen::Matrix4d::Identity();
    A(0, 1) = dt;
    A(2, 3) = dt;

    Eigen::Matrix<double, 4, 2> B = Eigen::Matrix<double, 4, 2>::Zero();
    B(0, 0) = 0.5 * dt * dt;
    B(1, 0) = dt;
    B(2, 1) = 0.5 * dt * dt;
    B(3, 1) = dt;

    Eigen::Matrix4d Q = Eigen::Matrix4d::Identity();
    Eigen::Matrix2d R = 0.1 * Eigen::Matrix2d::Identity();

    // ---- ctrlpp warm-up ----
    auto K_ctrlpp = ctrlpp::lqr_gain<double, NX, NU>(A, B, Q, R);

    // ---- Drake warm-up ----
    auto drake_result = drake::systems::controllers::LinearQuadraticRegulator(A, B, Q, R);

    // ---- Benchmark ----
    ankerl::nanobench::Bench bench;
    bench.title("LQR: ctrlpp vs Drake")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::lqr_gain",
             [&]
             {
                 auto K = ctrlpp::lqr_gain<double, NX, NU>(A, B, Q, R);
                 ankerl::nanobench::doNotOptimizeAway(K);
             })
        .run("drake::LinearQuadraticRegulator",
             [&]
             {
                 auto r = drake::systems::controllers::LinearQuadraticRegulator(A, B, Q, R);
                 ankerl::nanobench::doNotOptimizeAway(r);
             });

    std::ofstream csv("bench_lqr_vs_drake.csv");
    bench.render(comma_csv_tpl, csv);
}
