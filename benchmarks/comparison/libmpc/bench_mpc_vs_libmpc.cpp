// Competitive benchmark: ctrlpp::mpc vs libmpc++ LMPC
// Problem: Linear MPC, double-integrator (NX=4, NU=2), horizon N=10, box constraints

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/state_space.h"

#include <cassert>
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
    // Template: <NX, NU, Ndu, NY, Npred, Nctrl>
    constexpr std::size_t NY = NX;
    constexpr std::size_t Ndu = 0;
    constexpr int Nctrl = N;
    mpc::LMPC<NX, NU, Ndu, NY, N, Nctrl> libmpc_ctrl;

    mpc::mat<NX, NX> A_mpc = A;
    mpc::mat<NX, NU> B_mpc = B;
    mpc::mat<NY, NX> C_mpc = mpc::mat<NY, NX>::Identity();
    libmpc_ctrl.setStateSpaceModel(A_mpc, B_mpc, C_mpc);

    // Weights: full matrices (NY x N) and (NU x Nctrl)
    mpc::mat<NY, N> O_weight;
    O_weight.colwise() = mpc::cvec<NY>::Ones();
    mpc::mat<NU, Nctrl> U_weight;
    U_weight.colwise() = 0.1 * mpc::cvec<NU>::Ones();
    mpc::mat<NU, Nctrl> DU_weight = mpc::mat<NU, Nctrl>::Zero();
    libmpc_ctrl.setObjectiveWeights(O_weight, U_weight, DU_weight);

    // Bounds: full matrices (NX x N) and (NU x Nctrl)
    mpc::mat<NX, N> X_min_mpc;
    X_min_mpc.colwise() = mpc::cvec<NX>::Constant(-5.0);
    mpc::mat<NX, N> X_max_mpc;
    X_max_mpc.colwise() = mpc::cvec<NX>::Constant(5.0);
    libmpc_ctrl.setStateBounds(X_min_mpc, X_max_mpc);

    mpc::mat<NU, Nctrl> U_min_mpc;
    U_min_mpc.colwise() = mpc::cvec<NU>::Constant(-1.0);
    mpc::mat<NU, Nctrl> U_max_mpc;
    U_max_mpc.colwise() = mpc::cvec<NU>::Constant(1.0);
    libmpc_ctrl.setInputBounds(U_min_mpc, U_max_mpc);

    mpc::cvec<NX> x0_mpc = x0;
    mpc::cvec<NY> yref_mpc = mpc::cvec<NY>::Zero();
    // Warm up
    libmpc_ctrl.optimize(x0_mpc, mpc::cvec<NU>::Zero());

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
                 auto r = libmpc_ctrl.optimize(x0_mpc, mpc::cvec<NU>::Zero());
                 ankerl::nanobench::doNotOptimizeAway(r);
             });

    std::ofstream csv("bench_mpc_vs_libmpc.csv");
    bench.render(comma_csv_tpl, csv);
}
