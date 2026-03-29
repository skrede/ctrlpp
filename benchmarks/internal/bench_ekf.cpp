#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/ekf.h"

#include <cmath>
#include <fstream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

namespace
{

// 4-state nonlinear dynamics: constant-turn-rate model
struct dynamics
{
    auto operator()(const ctrlpp::Vector<double, 4>& x,
                    const ctrlpp::Vector<double, 1>& /*u*/) const -> ctrlpp::Vector<double, 4>
    {
        constexpr double dt = 0.01;
        constexpr double omega = 0.1;
        ctrlpp::Vector<double, 4> xn;
        xn[0] = x[0] + x[1] * dt;
        xn[1] = x[1] + std::cos(omega * dt) * 0.01;
        xn[2] = x[2] + x[3] * dt;
        xn[3] = x[3] + std::sin(omega * dt) * 0.01;
        return xn;
    }
};

// Observe position components
struct measurement
{
    auto operator()(const ctrlpp::Vector<double, 4>& x) const -> ctrlpp::Vector<double, 2>
    {
        return ctrlpp::Vector<double, 2>{x[0], x[2]};
    }
};

}

int main()
{
    ctrlpp::ekf_config<double, 4, 1, 2> cfg{};
    cfg.Q = Eigen::Matrix4d::Identity() * 0.01;
    cfg.R = Eigen::Matrix2d::Identity() * 0.1;
    cfg.x0 << 0.0, 1.0, 0.0, 1.0;

    ctrlpp::ekf filter(dynamics{}, measurement{}, cfg);

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    ctrlpp::Vector<double, 2> z{1.0, 0.5};

    ankerl::nanobench::Bench bench;
    bench.title("EKF")
        .warmup(50)
        .minEpochIterations(1000)
        .performanceCounters(true);

    bench.run("ekf::predict", [&] {
        filter.predict(u);
        ankerl::nanobench::doNotOptimizeAway(filter.state());
    });

    bench.run("ekf::update", [&] {
        filter.update(z);
        ankerl::nanobench::doNotOptimizeAway(filter.state());
    });

    std::ofstream csv("bench_ekf.csv");
    bench.render(csv_tpl, csv);
}
