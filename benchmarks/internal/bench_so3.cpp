#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "ctrlpp/lie/so3.h"

#include <cmath>
#include <fstream>

namespace
{

// Antipodal coefficient vectors name the same rotation, so the distance between
// two unit quaternions is the smaller of the two alignments.
double quaternion_deviation(const Eigen::Quaterniond& lhs, const Eigen::Quaterniond& rhs)
{
    return std::min((lhs.coeffs() - rhs.coeffs()).norm(), (lhs.coeffs() + rhs.coeffs()).norm());
}

}

int main(int argc, char** argv)
{
    // Rotation vector: 45 degrees about [1,1,1]/sqrt(3)
    ctrlpp::Vector<double, 3> omega;
    omega << 0.4534498, 0.4534498, 0.4534498;

    auto q = ctrlpp::so3::exp(omega);

    // Along one axis the exponential is a group homomorphism, so exp(2w) must
    // equal exp(w) composed with itself. Neither map encodes that identity: the
    // left side is one closed-form evaluation and the right is a Hamilton
    // product of two others, so what it leaves is the arithmetic.
    const double exp_residual = quaternion_deviation(
        ctrlpp::so3::exp(ctrlpp::Vector<double, 3>{2.0 * omega}), ctrlpp::so3::compose(q, q));
    const double log_residual = (ctrlpp::so3::log(q) - omega).norm() / omega.norm();

    ankerl::nanobench::Bench bench;
    bench.title("SO3")
        .warmup(100)
        .minEpochIterations(10000)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_single_implementation_row(bench, "so3::exp", [&] {
        auto result = ctrlpp::so3::exp(omega);
        ankerl::nanobench::doNotOptimizeAway(result);
    });

    ctrlpp::bench::run_single_implementation_row(bench, "so3::log", [&] {
        auto result = ctrlpp::so3::log(q);
        ankerl::nanobench::doNotOptimizeAway(result);
    });

    ctrlpp::bench::run_own_criterion_row(
        bench, "quaternion distance between this map's exp(2w) and its exp(w) composed with itself",
        "so3::exp", exp_residual, [&] {
            auto result = ctrlpp::so3::exp(omega);
            ankerl::nanobench::doNotOptimizeAway(result);
        });

    ctrlpp::bench::run_own_criterion_row(
        bench, "relative deviation of this map's log(exp(w)) from the rotation vector w it was built from",
        "so3::log", log_residual, [&] {
            auto result = ctrlpp::so3::log(q);
            ankerl::nanobench::doNotOptimizeAway(result);
        });

    std::ofstream csv("bench_so3.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
