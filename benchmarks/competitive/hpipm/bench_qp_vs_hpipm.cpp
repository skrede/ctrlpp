// Competitive benchmark: ctrlpp::mpc vs HPIPM (via hpipm-cpp wrapper)
// Problem: OCP-structured QP, double-integrator (NX=4, NU=2), horizon N=10
//
// Note: hpipm-cpp wrapper adds overhead compared to raw HPIPM C API.
// This comparison measures the full user-facing API path for both libraries.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/state_space.h"

#include <hpipm-cpp/hpipm-cpp.hpp>

#include <Eigen/Dense>

#include <fstream>
#include <vector>

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
    ctrlpp_mpc.solve(x0);

    // ---- HPIPM setup via hpipm-cpp ----
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    std::vector<int> nx_vec(N + 1, nx);
    std::vector<int> nu_vec(N + 1, nu);
    nu_vec[N] = 0; // no input at terminal stage

    hpipm::OcpQpDim dim;
    dim.N = N;
    dim.nx = nx_vec;
    dim.nu = nu_vec;

    // Set up OCP QP data
    hpipm::OcpQp qp(dim);
    for (int k = 0; k < N; ++k) {
        qp.A[k] = A;
        qp.B[k] = B;
        qp.Q[k] = Q;
        qp.R[k] = R;
        qp.q[k] = Eigen::Vector4d::Zero();
        qp.r[k] = Eigen::Vector2d::Zero();

        qp.idxbu[k] = {0, 1};
        qp.lbu[k] = Eigen::Vector2d::Constant(-1.0);
        qp.ubu[k] = Eigen::Vector2d::Constant(1.0);

        qp.idxbx[k] = {0, 1, 2, 3};
        qp.lbx[k] = Eigen::Vector4d::Constant(-5.0);
        qp.ubx[k] = Eigen::Vector4d::Constant(5.0);
    }
    qp.Q[N] = Q;
    qp.q[N] = Eigen::Vector4d::Zero();

    hpipm::OcpQpIpmSolverSettings settings;
    settings.mode = hpipm::HpipmMode::Balance;
    settings.iter_max = 30;
    settings.tol_stat = 1e-6;
    settings.tol_eq = 1e-6;
    settings.tol_ineq = 1e-6;
    settings.tol_comp = 1e-6;

    hpipm::OcpQpSolver solver(dim, settings);
    hpipm::OcpQpSolution solution(dim);

    // Warm up
    solver.solve(x0, qp, solution);

    // ---- Benchmark ----
    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp vs HPIPM")
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
        .run("hpipm::OcpQpSolver::solve",
             [&]
             {
                 solver.solve(x0, qp, solution);
                 auto u0 = solution.u[0];
                 ankerl::nanobench::doNotOptimizeAway(u0);
             });

    std::ofstream csv("bench_qp_vs_hpipm.csv");
    bench.render(comma_csv_tpl, csv);
}
