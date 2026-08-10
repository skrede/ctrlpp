#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/kalman.h"
#include "ctrlpp/model/state_space.h"

#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
    // 4-state, 2-measurement, 1-input system
    constexpr double dt = 0.01;

    ctrlpp::discrete_state_space<double, 4, 1, 2> sys{};
    sys.A << 1.0, dt, 0.0, 0.0,
             0.0, 1.0, dt, 0.0,
             0.0, 0.0, 1.0, dt,
             0.0, 0.0, 0.0, 1.0;
    sys.B << 0.0, dt, 0.0, 0.0;
    sys.C << 1.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 0.0;
    sys.D.setZero();

    ctrlpp::kalman_config<double, 4, 1, 2> cfg{};
    cfg.Q = Eigen::Matrix4d::Identity() * 0.01;
    cfg.R = Eigen::Matrix2d::Identity() * 0.1;

    auto kf = ctrlpp::bench::built_or_exit(ctrlpp::kalman_filter<double, 4, 1, 2>::create(sys, cfg), "kalman_filter");

    Eigen::Matrix<double, 1, 1> u;
    u << 0.5;
    Eigen::Vector2d z{1.0, 0.5};

    ankerl::nanobench::Bench bench;
    bench.title("Kalman")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_single_implementation_row(bench, "kf::predict", [&] {
        kf.predict(u);
        ankerl::nanobench::doNotOptimizeAway(kf.state());
    });

    // A rejected step performs a fraction of the work, so a run that
    // included one would report a meaningless figure. Establish outside the
    // measured region that the step runs, then feed the result to the
    // optimizer barrier inside it so it is neither discarded nor branched on
    // while the clock is running.
    if(const auto stepped = kf.update(z); !stepped)
    {
        std::cerr << "bench_kalman: kalman_filter::update rejected the measurement; the reported figures would be meaningless\n";
        return 1;
    }

    ctrlpp::bench::run_single_implementation_row(bench, "kf::update", [&] {
        ankerl::nanobench::doNotOptimizeAway(kf.update(z).has_value());
        ankerl::nanobench::doNotOptimizeAway(kf.state());
    });

    std::ofstream csv("bench_kalman.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
