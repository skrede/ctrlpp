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

#include "ctrlpp/control/care.h"
#include "ctrlpp/control/lqr.h"

#include "bench_construct.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/core/types/StateVector.h>
#include <ct/core/types/ControlVector.h>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>
#include <ct/optcon/lqr/LQR.hpp>
#include <ct/optcon/lqr/LQR-impl.hpp>

#include <Eigen/Dense>

#include <fstream>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

template <std::size_t NX, std::size_t NU>
auto build_damped_chain()
{
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Zero();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) = -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = 1.0;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = 1.0;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();

    return std::tuple{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    auto [A, B, Q, R] = build_damped_chain<NX, NU>();

    // The warmup also asserts the solve succeeds: a benchmark that times a
    // refused solve reports a number for a problem the library declined.
    (void)ctrlpp::bench::built_or_exit(ctrlpp::lqr_gain_continuous<double, NX, NU>(A, B, Q, R),
                                       "lqr_gain_continuous warmup on the damped chain");

    ct::optcon::LQR<NX, NU> ct_lqr;
    typename ct::optcon::LQR<NX, NU>::state_matrix_t A_ct = A;
    typename ct::optcon::LQR<NX, NU>::state_matrix_t Q_ct = Q;
    typename ct::optcon::LQR<NX, NU>::control_matrix_t R_ct = R;
    Eigen::Matrix<double, int(NX), int(NU)> B_ct = B;
    Eigen::Matrix<double, int(NU), int(NX)> K_ct;

    ct_lqr.compute(Q_ct, R_ct, A_ct, B_ct, K_ct);

    bench.run(label_ctrlpp,
              [&]
              {
                  auto K = ctrlpp::lqr_gain_continuous<double, NX, NU>(A, B, Q, R);
                  ankerl::nanobench::doNotOptimizeAway(K);
              })
        .run(label_ct,
             [&]
             {
                 ct_lqr.compute(Q_ct, R_ct, A_ct, B_ct, K_ct);
                 ankerl::nanobench::doNotOptimizeAway(K_ct);
             });
}

} // namespace

int main()
{
    ankerl::nanobench::Bench bench;
    bench.title("Continuous LQR: ctrlpp::lqr_gain_continuous vs ct::optcon::LQR (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);

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
    bench.render(comma_csv_tpl, csv);
}
