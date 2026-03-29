// Competitive benchmark: ctrlpp::mpc vs libmpc++ LMPC
// Problem: Linear MPC, double-integrator (NX=4, NU=2), horizon N=10, box constraints

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/state_space.h"

#include <mpc/LMPC.hpp>

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
    // ---- Problem setup: discrete double-integrator, NX=4, NU=2 ----
    // x = [p1, v1, p2, v2], u = [a1, a2]
    // x_{k+1} = A x_k + B u_k, dt = 0.1
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;
    constexpr int N = 10;
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

    Eigen::Vector4d x0;
    x0 << 1.0, 0.0, -0.5, 0.0;

    // ---- ctrlpp MPC setup ----
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = A,
        .B = B,
        .C = Eigen::Matrix4d::Identity(),
        .D = Eigen::Matrix<double, 4, 2>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{};
    cfg.horizon = N;
    cfg.Q = Q;
    cfg.R = R;
    cfg.u_min = Eigen::Vector2d::Constant(-1.0);
    cfg.u_max = Eigen::Vector2d::Constant(1.0);
    cfg.x_min = Eigen::Vector4d::Constant(-5.0);
    cfg.x_max = Eigen::Vector4d::Constant(5.0);

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> ctrlpp_mpc(sys, cfg);
    // Warm up solver
    ctrlpp_mpc.solve(x0);

    // ---- libmpc++ setup ----
    mpc::LMPC<NX, NU, 0, N> libmpc_ctrl;

    mpc::mat<NX, NX> A_mpc = A;
    mpc::mat<NX, NU> B_mpc = B;
    mpc::mat<NX, NX> C_mpc = mpc::mat<NX, NX>::Identity();
    libmpc_ctrl.setStateSpaceModel(A_mpc, B_mpc, C_mpc);

    mpc::mat<NX, NX> Q_mpc = Q;
    mpc::mat<NU, NU> R_mpc = R;
    libmpc_ctrl.setObjectiveWeights(Q_mpc, R_mpc, Q_mpc);

    mpc::cvec<NU> u_min_mpc;
    u_min_mpc << -1.0, -1.0;
    mpc::cvec<NU> u_max_mpc;
    u_max_mpc << 1.0, 1.0;
    libmpc_ctrl.setInputBounds(u_min_mpc, u_max_mpc);

    mpc::cvec<NX> x_min_mpc = mpc::cvec<NX>::Constant(-5.0);
    mpc::cvec<NX> x_max_mpc = mpc::cvec<NX>::Constant(5.0);
    libmpc_ctrl.setStateBounds(x_min_mpc, x_max_mpc);

    mpc::cvec<NX> x0_mpc = x0;
    // Warm up
    libmpc_ctrl.step(x0_mpc, mpc::cvec<NX>::Zero());

    // ---- Benchmark ----
    ankerl::nanobench::Bench bench;
    bench.title("MPC: ctrlpp vs libmpc++")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::mpc::solve",
             [&]
             {
                 auto u = ctrlpp_mpc.solve(x0);
                 ankerl::nanobench::doNotOptimizeAway(u);
             })
        .run("libmpc++::LMPC::step",
             [&]
             {
                 auto r = libmpc_ctrl.step(x0_mpc, mpc::cvec<NX>::Zero());
                 ankerl::nanobench::doNotOptimizeAway(r);
             });

    std::ofstream csv("bench_mpc_vs_libmpc.csv");
    bench.render(comma_csv_tpl, csv);
}
