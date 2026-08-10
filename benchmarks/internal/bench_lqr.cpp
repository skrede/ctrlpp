#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "comparison/ct/riccati_problem.h"

#include "ctrlpp/control/lqr.h"

#include <fstream>

int main(int argc, char** argv)
{
    // Double integrator: x1(k+1) = x1(k) + dt*x2(k), x2(k+1) = x2(k) + dt*u(k)
    constexpr double dt = 0.01;
    ctrlpp::bench::riccati_plant<2, 1> plant{};
    plant.A << 1.0, dt, 0.0, 1.0;
    plant.B << 0.0, dt;
    plant.Q = Eigen::Matrix2d::Identity();
    plant.R = Eigen::Matrix<double, 1, 1>::Identity();

    const Eigen::Matrix<double, 1, 2> K = ctrlpp::bench::built_or_exit(
        ctrlpp::lqr_gain<double, 2, 1>(plant.A, plant.B, plant.Q, plant.R),
        "lqr_gain on the double integrator");
    ctrlpp::lqr<double, 2, 1> controller(K);
    Eigen::Vector2d x{1.0, 0.5};

    ankerl::nanobench::Bench bench;
    bench.title("LQR")
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    bench.warmup(50).minEpochIterations(1000);
    ctrlpp::bench::run_single_implementation_row(bench, "lqr_gain", [&] {
        auto gain = ctrlpp::lqr_gain<double, 2, 1>(plant.A, plant.B, plant.Q, plant.R);
        ankerl::nanobench::doNotOptimizeAway(gain);
    });

    bench.warmup(100).minEpochIterations(10000);
    ctrlpp::bench::run_single_implementation_row(bench, "lqr::compute", [&] {
        auto u = controller.compute(x);
        ankerl::nanobench::doNotOptimizeAway(u);
    });

    bench.warmup(50).minEpochIterations(1000);
    ctrlpp::bench::run_own_criterion_row(
        bench, "relative residual of the discrete Riccati solution reconstructed from this arm's own gain",
        "lqr_gain", ctrlpp::bench::discrete_gain_optimality_residual<2, 1>(plant, K), [&] {
            auto gain = ctrlpp::lqr_gain<double, 2, 1>(plant.A, plant.B, plant.Q, plant.R);
            ankerl::nanobench::doNotOptimizeAway(gain);
        });

    ctrlpp::bench::run_certificate_row(
        bench, "spectral radius of the discrete closed loop this arm's gain forms (below one certifies stability)",
        "lqr_gain", ctrlpp::bench::closed_loop_spectral_radius<2, 1>(plant, K), [&] {
            auto gain = ctrlpp::lqr_gain<double, 2, 1>(plant.A, plant.B, plant.Q, plant.R);
            ankerl::nanobench::doNotOptimizeAway(gain);
        });

    std::ofstream csv("bench_lqr.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
