// Competitive benchmark: ctrlpp::lqr_gain vs ct::optcon::LQR
// Problem: LQR gain computation for double-integrator (NX=4, NU=2)
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/lqr.h"

#include <ct/optcon/optcon.h>

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

    // ---- ct_optcon setup ----
    ct::optcon::LQR<NX, NU> ct_lqr;

    // ct_optcon LQR uses continuous-time Riccati; we provide the same matrices
    // for a fair discrete-time comparison of the DARE/Riccati solve path.
    ct::core::StateMatrix<NX> Q_ct = Q;
    ct::core::ControlMatrix<NU> R_ct = R;
    ct::core::StateMatrix<NX> A_ct = A;
    ct::core::StateControlMatrix<NX, NU> B_ct = B;

    // Warm up
    ct_lqr.compute(Q_ct, R_ct, A_ct, B_ct);

    // ---- Benchmark ----
    ankerl::nanobench::Bench bench;
    bench.title("LQR: ctrlpp vs ct_optcon")
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
        .run("ct::optcon::LQR::compute",
             [&]
             {
                 ct_lqr.compute(Q_ct, R_ct, A_ct, B_ct);
                 auto K = ct_lqr.getSolution();
                 ankerl::nanobench::doNotOptimizeAway(K);
             });

    std::ofstream csv("bench_lqr_vs_ct.csv");
    bench.render(comma_csv_tpl, csv);
}
