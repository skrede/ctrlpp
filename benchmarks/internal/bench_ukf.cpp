#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/ukf.h"

#include <cmath>
#include <fstream>
#include <iostream>

static constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

namespace
{

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
    ctrlpp::ukf_config<double, 4, 1, 2> cfg{};
    cfg.Q = Eigen::Matrix4d::Identity() * 0.01;
    cfg.R = Eigen::Matrix2d::Identity() * 0.1;
    cfg.x0 << 0.0, 1.0, 0.0, 1.0;

    ctrlpp::ukf filter(dynamics{}, measurement{}, cfg);

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    ctrlpp::Vector<double, 2> z{1.0, 0.5};

    ankerl::nanobench::Bench bench;
    bench.title("UKF")
        .warmup(50)
        .minEpochIterations(1000)
        .performanceCounters(true);

    bench.run("ukf::predict", [&] {
        filter.predict(u);
        ankerl::nanobench::doNotOptimizeAway(filter.state());
    });

    // A rejected step performs a fraction of the work, so a run that
    // included one would report a meaningless figure. Establish outside the
    // measured region that the step runs, then feed the result to the
    // optimizer barrier inside it so it is neither discarded nor branched on
    // while the clock is running.
    if(const auto stepped = filter.update(z); !stepped)
    {
        std::cerr << "bench_ukf: ukf::update rejected the measurement; the reported figures would be meaningless\n";
        return 1;
    }

    bench.run("ukf::update", [&] {
        ankerl::nanobench::doNotOptimizeAway(filter.update(z).has_value());
        ankerl::nanobench::doNotOptimizeAway(filter.state());
    });

    std::ofstream csv("bench_ukf.csv");
    bench.render(csv_tpl, csv);
}
