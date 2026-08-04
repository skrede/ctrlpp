// Discrete Riccati benchmark, swept over state dimension.
//
// The acceptance check is the independent variable. Three quantities are timed
// per size on identical operands:
//
//   solve+accept -- the whole `ctrlpp::dare` entry point
//   accept       -- `verify_dare_solution` alone, on the solve's own answer
//   fwd-error    -- `estimate_riccati_forward_error` alone, on pre-formed
//                   closed-loop and residual operands
//
// The solve without the check is the difference of the first two, which is the
// denominator the published operation counts use: those are stated against a
// solve with no acceptance arithmetic in it. Timing the round's start against
// its end would confound the construction change, the removed factorization and
// the operand preparation into one number that attributes to none of them.
//
// The corpus is the discrete damped chain: a forward-Euler step of the same
// continuous chain the CARE bakeoff sweeps, so the two harnesses describe the
// same family. Q = I and R = 0.1 I put the weight scale at exactly one, so
// equilibration is inactive and the solve runs the acceptance check ONCE. On a
// pose whose weights are not already equilibrated the entry point checks twice,
// at the equilibrated scale and again at the caller's, plus a definiteness
// factorization and a gain agreement; that path is not what these rows measure.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/control/dare.h"

#include <Eigen/Dense>

#include <cstdio>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <tuple>
#include <algorithm>

namespace
{

constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

/// Forward-Euler step of the continuous damped chain the CARE bakeoff uses:
/// -0.5 on the diagonal and 1.0 on the superdiagonal, one input per group of
/// states. The step keeps every eigenvalue of A at 1 - 0.5 dt, inside the unit
/// disk, so the discrete equation is well posed at every size.
template <std::size_t NX, std::size_t NU>
auto build_discrete_damped_chain()
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);
    constexpr double dt = 0.01;

    Eigen::Matrix<double, n, n> A = Eigen::Matrix<double, n, n>::Identity();
    for (std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5;
    for (std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt * 1.0;

    Eigen::Matrix<double, n, nu> B = Eigen::Matrix<double, n, nu>::Zero();
    const std::size_t group = NX / NU;
    for (std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = dt;
    }

    Eigen::Matrix<double, n, n>   Q = Eigen::Matrix<double, n, n>::Identity();
    Eigen::Matrix<double, nu, nu> R = 0.1 * Eigen::Matrix<double, nu, nu>::Identity();

    return std::tuple{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench)
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);

    auto [A, B, Q, R] = build_discrete_damped_chain<NX, NU>();

    // Solved once, outside every measurement window, so that the acceptance
    // rows are handed the same answer the solve produced rather than one they
    // recompute.
    auto solved = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
    if (!solved)
    {
        std::cerr << "SKIP NX=" << NX << ": dare declined the corpus\n";
        return;
    }
    const Eigen::Matrix<double, n, n>  P = solved->P;
    const Eigen::Matrix<double, nu, n> K = solved->K;

    // The two operands the forward-error estimate consumes, formed exactly as
    // the acceptance path forms them. Building them here rather than inside the
    // timed loop is what makes the fwd-error row the estimate's own cost and
    // not the cost of preparing its inputs.
    const Eigen::Matrix<double, n, n> AtPA        = (A.transpose() * P * A).eval();
    const Eigen::Matrix<double, n, n> AtPBK       = (A.transpose() * P * B * K).eval();
    const Eigen::Matrix<double, n, n> residual    = (AtPA - P - AtPBK + Q).eval();
    const Eigen::Matrix<double, n, n> closed_loop = (A - B * K).eval();

    char buf[64];

    std::snprintf(buf, sizeof(buf), "dare solve+accept NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto r = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
        ankerl::nanobench::doNotOptimizeAway(r);
    });

    std::snprintf(buf, sizeof(buf), "dare accept NX=%zu", NX);
    bench.run(buf, [&]
    {
        Eigen::Matrix<double, nu, n> verified_gain;
        auto v = ctrlpp::detail::verify_dare_solution<double, NX, NU>(A, B, Q, R, P, &verified_gain);
        ankerl::nanobench::doNotOptimizeAway(v);
        ankerl::nanobench::doNotOptimizeAway(verified_gain);
    });

    std::snprintf(buf, sizeof(buf), "dare fwd-error NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto e = ctrlpp::detail::estimate_riccati_forward_error<double, n>(closed_loop, residual, P);
        ankerl::nanobench::doNotOptimizeAway(e);
    });
}

void check_perf_event_paranoid()
{
    std::ifstream paranoid_file("/proc/sys/kernel/perf_event_paranoid");
    int paranoid_value = 3;
    paranoid_file >> paranoid_value;
    if (paranoid_value > 2)
        std::cerr << "WARN: perf_event_paranoid=" << paranoid_value
                  << "; nanobench hardware counters may be zero. "
                  << "Run:  sudo sysctl kernel.perf_event_paranoid=2\n";
}

}

int main()
{
    check_perf_event_paranoid();

    // 51 epochs give an odd-sized sample, so the reported median is an observed
    // measurement rather than an interpolation, and the CSV's error column is
    // the median absolute percent error across those same 51.
    ankerl::nanobench::Bench bench;
    bench.title("DARE acceptance overhead (size sweep)")
        .warmup(10)
        .epochs(51)
        .minEpochTime(std::chrono::milliseconds(1))
        .performanceCounters(true);

    // NX = 15 IS THE CEILING, AND IT IS A COMPILE-TIME ONE. The forward-error
    // estimate holds an M-by-M operator with M = NX(NX+1)/2 as a fixed-size
    // Eigen object, so Eigen's 128 KiB stack-allocation limit admits M <= 128,
    // that is NX <= 15. At NX = 16 the operator is 136 x 136 = 147,968 bytes and
    // the static assertion fires: the acceptance check cannot be instantiated
    // there at all, so there is no larger size for this sweep to reach.
    run_size_sweep<2,  1>(bench);
    run_size_sweep<4,  2>(bench);
    run_size_sweep<6,  2>(bench);
    run_size_sweep<8,  2>(bench);
    run_size_sweep<12, 3>(bench);
    run_size_sweep<15, 3>(bench);

    std::ofstream csv("bench_dare.csv");
    bench.render(csv_tpl, csv);
}
