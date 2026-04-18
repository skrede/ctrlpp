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

} // namespace

int main()
{
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;
    constexpr double dt = 0.1;

    // Same problem setup as bench_lqr_vs_ct.cpp for a fair head-to-head.
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

    // ---- ctrlpp warm-up: continuous-time LQR via care + K = R^-1 B^T P ----
    auto warmup_care = ctrlpp::care<double, NX, NU>(A, B, Q, R);
    // Accept failure on warm-up; the benchmark will still run the timing loop.
    (void)warmup_care;

    // ---- ct_optcon setup ----
    ct::optcon::LQR<NX, NU> ct_lqr;
    using LQR_t = ct::optcon::LQR<NX, NU>;
    typename LQR_t::state_matrix_t Q_ct = Q;
    typename LQR_t::control_matrix_t R_ct = R;
    typename LQR_t::state_matrix_t A_ct = A;
    Eigen::Matrix<double, NX, NU> B_ct = B;
    Eigen::Matrix<double, NU, NX> K_ct;

    ct_lqr.compute(Q_ct, R_ct, A_ct, B_ct, K_ct);

    // ---- Benchmark: apples-to-apples continuous-time LQR ----
    ankerl::nanobench::Bench bench;
    bench.title("Continuous LQR: ctrlpp (care + K) vs ct_optcon")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::lqr_gain_continuous",
             [&]
             {
                 auto K = ctrlpp::lqr_gain_continuous<double, NX, NU>(A, B, Q, R);
                 ankerl::nanobench::doNotOptimizeAway(K);
             })
        .run("ct::optcon::LQR::compute",
             [&]
             {
                 ct_lqr.compute(Q_ct, R_ct, A_ct, B_ct, K_ct);
                 ankerl::nanobench::doNotOptimizeAway(K_ct);
             });

    std::ofstream csv("bench_lqr_continuous_vs_ct.csv");
    bench.render(comma_csv_tpl, csv);
}
