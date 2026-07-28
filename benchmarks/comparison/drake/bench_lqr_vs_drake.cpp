// Competitive benchmark: ctrlpp::lqr_gain (DARE-based) vs Drake LinearQuadraticRegulator
// Size-swept NX in {2, 4, 6, 8, 12, 16, 20, 24, 30}. NU scales with NX.
//
// Drake requires manual installation. No AUR package. Bazel-only build system.
// Pre-built tar.gz available for Ubuntu/macOS. See benchmarks/README.md.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/lqr.h"

#include "bench_construct.h"

#include <drake/systems/controllers/linear_quadratic_regulator.h>

#include <Eigen/Dense>

#include <fstream>
#include <cstddef>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

template <std::size_t NX, std::size_t NU>
auto build_chain_of_integrators(double dt)
{
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = dt;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();

    return std::tuple{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_drake)
{
    constexpr double dt = 0.05;
    auto [A, B, Q, R] = build_chain_of_integrators<NX, NU>(dt);

    // The warmup also asserts the solve succeeds: a benchmark that times a
    // refused solve reports a number for a problem the library declined.
    (void)ctrlpp::bench::built_or_exit(ctrlpp::lqr_gain<double, NX, NU>(A, B, Q, R),
                                       "lqr_gain warmup on the chain of integrators");

    auto warmup_drake = drake::systems::controllers::LinearQuadraticRegulator(A, B, Q, R);
    (void)warmup_drake;

    bench.run(label_ctrlpp,
              [&]
              {
                  auto K = ctrlpp::lqr_gain<double, NX, NU>(A, B, Q, R);
                  ankerl::nanobench::doNotOptimizeAway(K);
              })
        .run(label_drake,
             [&]
             {
                 auto r = drake::systems::controllers::LinearQuadraticRegulator(A, B, Q, R);
                 ankerl::nanobench::doNotOptimizeAway(r);
             });
}

} // namespace

int main()
{
    ankerl::nanobench::Bench bench;
    bench.title("LQR: ctrlpp::lqr_gain (DARE) vs drake::LinearQuadraticRegulator (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);

    run_size_sweep<2, 1>(bench,  "ctrlpp::lqr_gain NX=2",  "drake::LQR NX=2");
    run_size_sweep<4, 2>(bench,  "ctrlpp::lqr_gain NX=4",  "drake::LQR NX=4");
    run_size_sweep<6, 2>(bench,  "ctrlpp::lqr_gain NX=6",  "drake::LQR NX=6");
    run_size_sweep<8, 2>(bench,  "ctrlpp::lqr_gain NX=8",  "drake::LQR NX=8");
    run_size_sweep<12, 3>(bench, "ctrlpp::lqr_gain NX=12", "drake::LQR NX=12");
    run_size_sweep<16, 4>(bench, "ctrlpp::lqr_gain NX=16", "drake::LQR NX=16");
    run_size_sweep<20, 5>(bench, "ctrlpp::lqr_gain NX=20", "drake::LQR NX=20");
    run_size_sweep<24, 6>(bench, "ctrlpp::lqr_gain NX=24", "drake::LQR NX=24");
    run_size_sweep<30, 6>(bench, "ctrlpp::lqr_gain NX=30", "drake::LQR NX=30");

    std::ofstream csv("bench_lqr_vs_drake.csv");
    bench.render(comma_csv_tpl, csv);
}
