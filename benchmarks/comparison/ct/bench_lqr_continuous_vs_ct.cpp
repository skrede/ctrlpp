// Competitive benchmark: ctrlpp continuous-time LQR (via care + K = R^-1 B^T P) vs ct::optcon::LQR
// Problem: identical to bench_lqr_vs_ct.cpp (integrator dynamics, NX=4, NU=2, dt=0.1)
//
// ct::optcon::LQR::compute is a continuous-time LQR built on ct::optcon::CARE::solve
// (Eigen::RealSchur + LAPACK dtrsen_ reorder). bench_lqr_vs_ct compares against
// ctrlpp::lqr_gain which chains through DARE -- apples to oranges. This benchmark
// chains ctrlpp::care for the continuous-time Riccati and computes K directly, giving
// a straight apples-to-apples race.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ct_lqr_arm.h"
#include "riccati_problem.h"

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/care.h"

#include <Eigen/Dense>

#include <cstddef>
#include <fstream>

namespace
{

using ctrlpp::bench::build_damped_chain;
using ctrlpp::bench::closed_loop_abscissa;
using ctrlpp::bench::ct_lqr_arm;
using ctrlpp::bench::gain_optimality_residual;
using ctrlpp::bench::riccati_plant;

constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' gains K";
constexpr char const* residual_metric =
    "relative residual of the Riccati solution reconstructed from this arm's own gain";
constexpr char const* abscissa_metric =
    "closed-loop spectral abscissa of this arm's own gain (negative certifies stability)";

// The residual says the gain is right and the abscissa says it is admissible;
// neither substitutes for the other, so both are published. Each reconstructed
// cost matrix comes from that arm's own gain rather than from a second Riccati
// solve, which would measure a different object than the one being benchmarked.
template <std::size_t NX, std::size_t NU>
void emit_rows(ankerl::nanobench::Bench& bench, const riccati_plant<NX, NU>& plant, ct_lqr_arm<NX, NU>& ct_arm,
               const char* label_ctrlpp, const char* label_ct,
               const Eigen::Matrix<double, int(NU), int(NX)>& K_ctrlpp,
               const Eigen::Matrix<double, int(NU), int(NX)>& K_ct)
{
    auto solve_ctrlpp = [&]
    {
        auto K = ctrlpp::lqr_gain_continuous<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R);
        ankerl::nanobench::doNotOptimizeAway(K);
    };
    auto solve_ct = [&]
    {
        ct_arm.solve();
        ankerl::nanobench::doNotOptimizeAway(ct_arm.gain());
    };

    ctrlpp::bench::report_accuracy(bench, deviation_metric, (K_ctrlpp - K_ct).cwiseAbs().maxCoeff());
    bench.run(label_ctrlpp, solve_ctrlpp).run(label_ct, solve_ct);
    ctrlpp::bench::run_own_criterion_pair(bench, residual_metric, label_ctrlpp,
                                          gain_optimality_residual<NX, NU>(plant, K_ctrlpp), solve_ctrlpp, label_ct,
                                          gain_optimality_residual<NX, NU>(plant, K_ct), solve_ct);
    ctrlpp::bench::run_certificate_pair(bench, abscissa_metric, label_ctrlpp,
                                        closed_loop_abscissa<NX, NU>(plant, K_ctrlpp), solve_ctrlpp, label_ct,
                                        closed_loop_abscissa<NX, NU>(plant, K_ct), solve_ct);
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    const riccati_plant<NX, NU> plant = build_damped_chain<NX, NU>();
    const Eigen::Matrix<double, int(NU), int(NX)> K_ctrlpp = ctrlpp::bench::built_or_exit(
        ctrlpp::lqr_gain_continuous<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R), label_ctrlpp);
    ct_lqr_arm<NX, NU> ct_arm{plant};
    ct_arm.solve();
    const Eigen::Matrix<double, int(NU), int(NX)> K_ct = ct_arm.gain();
    emit_rows<NX, NU>(bench, plant, ct_arm, label_ctrlpp, label_ct, K_ctrlpp, K_ct);
}

} // namespace

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("Continuous LQR: ctrlpp::lqr_gain_continuous vs ct::optcon::LQR (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    run_size_sweep<2, 1>(bench,  "ctrlpp::lqr_gain_continuous NX=2",  "ct::optcon::LQR NX=2");
    run_size_sweep<4, 2>(bench,  "ctrlpp::lqr_gain_continuous NX=4",  "ct::optcon::LQR NX=4");
    run_size_sweep<6, 2>(bench,  "ctrlpp::lqr_gain_continuous NX=6",  "ct::optcon::LQR NX=6");
    run_size_sweep<8, 2>(bench,  "ctrlpp::lqr_gain_continuous NX=8",  "ct::optcon::LQR NX=8");
    run_size_sweep<12, 3>(bench, "ctrlpp::lqr_gain_continuous NX=12", "ct::optcon::LQR NX=12");
    run_size_sweep<16, 4>(bench, "ctrlpp::lqr_gain_continuous NX=16", "ct::optcon::LQR NX=16");
    run_size_sweep<20, 5>(bench, "ctrlpp::lqr_gain_continuous NX=20", "ct::optcon::LQR NX=20");
    run_size_sweep<24, 6>(bench, "ctrlpp::lqr_gain_continuous NX=24", "ct::optcon::LQR NX=24");
    run_size_sweep<30, 6>(bench, "ctrlpp::lqr_gain_continuous NX=30", "ct::optcon::LQR NX=30");

    std::ofstream csv("bench_lqr_continuous_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
