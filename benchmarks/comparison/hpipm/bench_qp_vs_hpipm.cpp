// Competitive benchmark: ctrlpp::mpc vs HPIPM (raw C API)
// Problem: OCP-structured QP, double-integrator (NX=4, NU=2), horizon N=10
//
// HPIPM is invoked via its C API directly (the upstream hpipm-cpp wrapper is
// out of sync with current HPIPM struct definitions). The thin adapter below
// allocates dim/qp/sol/arg/ws via the canonical memsize+create pattern.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_construct.h"

#include "lmpc/double_integrator.h"

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/state_space.h"

extern "C"
{
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_common.h>
}

#include <Eigen/Dense>

#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <memory>
#include <vector>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

struct aligned_buffer
{
    void* p = nullptr;

    explicit aligned_buffer(std::size_t bytes)
    {
        if(bytes == 0)
            return;
        p = std::aligned_alloc(64, ((bytes + 63u) / 64u) * 64u);
    }

    aligned_buffer(const aligned_buffer&) = delete;
    aligned_buffer& operator=(const aligned_buffer&) = delete;
    aligned_buffer(aligned_buffer&& o) noexcept : p{o.p} { o.p = nullptr; }
    aligned_buffer& operator=(aligned_buffer&& o) noexcept
    {
        if(this != &o) { std::free(p); p = o.p; o.p = nullptr; }
        return *this;
    }
    ~aligned_buffer() { std::free(p); }
};

}

int main()
{
    namespace problems = ctrlpp::bench::problems::lmpc;

    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;
    constexpr int N = 10;

    auto sys = problems::make_double_integrator_4_2_state_space();
    auto cfg = problems::make_double_integrator_4_2_config(N);
    auto x0 = problems::double_integrator_4_2_x0_default();

    auto ctrlpp_mpc = ctrlpp::bench::built_or_exit(
        ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, cfg), "ctrlpp_mpc");
    [[maybe_unused]] auto warm = ctrlpp_mpc.solve(x0);

    // --- HPIPM dimensions -----------------------------------------------------
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    std::vector<int> nx_vec(N + 1, nx);
    std::vector<int> nu_vec(N + 1, nu);
    nu_vec[N] = 0;
    std::vector<int> nbx_vec(N + 1, nx);
    std::vector<int> nbu_vec(N + 1, nu);
    nbu_vec[N] = 0;
    std::vector<int> ng_vec(N + 1, 0);
    std::vector<int> ns_vec(N + 1, 0);
    std::vector<int> nbxe_vec(N + 1, 0);
    nbxe_vec[0] = nx;

    aligned_buffer dim_buf{static_cast<std::size_t>(d_ocp_qp_dim_memsize(N))};
    d_ocp_qp_dim dim{};
    d_ocp_qp_dim_create(N, &dim, dim_buf.p);
    d_ocp_qp_dim_set_all(nx_vec.data(), nu_vec.data(),
                         nbx_vec.data(), nbu_vec.data(),
                         ng_vec.data(),  ns_vec.data(), &dim);
    for(int k = 0; k <= N; ++k)
        d_ocp_qp_dim_set_nbxe(k, nbxe_vec[k], &dim);

    // --- HPIPM QP data --------------------------------------------------------
    aligned_buffer qp_buf{static_cast<std::size_t>(d_ocp_qp_memsize(&dim))};
    d_ocp_qp qp{};
    d_ocp_qp_create(&dim, &qp, qp_buf.p);
    d_ocp_qp_set_all_zero(&qp);

    Eigen::Matrix<double, NX, NX> A_stage = sys.A;
    Eigen::Matrix<double, NX, NU> B_stage = sys.B;
    Eigen::Vector<double, NX> b_stage = Eigen::Vector<double, NX>::Zero();
    Eigen::Matrix<double, NX, NX> Q_stage = cfg.Q;
    Eigen::Matrix<double, NU, NU> R_stage = cfg.R;
    Eigen::Matrix<double, NU, NX> S_stage = Eigen::Matrix<double, NU, NX>::Zero();
    Eigen::Vector<double, NX> q_stage = Eigen::Vector<double, NX>::Zero();
    Eigen::Vector<double, NU> r_stage = Eigen::Vector<double, NU>::Zero();

    Eigen::Vector<double, NX> lbx_stage = *cfg.x_min;
    Eigen::Vector<double, NX> ubx_stage = *cfg.x_max;
    Eigen::Vector<double, NU> lbu_stage = *cfg.u_min;
    Eigen::Vector<double, NU> ubu_stage = *cfg.u_max;
    std::vector<int> idxbx{0, 1, 2, 3};
    std::vector<int> idxbu{0, 1};

    for(int k = 0; k < N; ++k)
    {
        d_ocp_qp_set_A(k, A_stage.data(), &qp);
        d_ocp_qp_set_B(k, B_stage.data(), &qp);
        d_ocp_qp_set_b(k, b_stage.data(), &qp);
        d_ocp_qp_set_Q(k, Q_stage.data(), &qp);
        d_ocp_qp_set_S(k, S_stage.data(), &qp);
        d_ocp_qp_set_R(k, R_stage.data(), &qp);
        d_ocp_qp_set_q(k, q_stage.data(), &qp);
        d_ocp_qp_set_r(k, r_stage.data(), &qp);

        d_ocp_qp_set_idxbx(k, idxbx.data(), &qp);
        d_ocp_qp_set_lbx(k, lbx_stage.data(), &qp);
        d_ocp_qp_set_ubx(k, ubx_stage.data(), &qp);
        d_ocp_qp_set_idxbu(k, idxbu.data(), &qp);
        d_ocp_qp_set_lbu(k, lbu_stage.data(), &qp);
        d_ocp_qp_set_ubu(k, ubu_stage.data(), &qp);
    }
    d_ocp_qp_set_Q(N, Q_stage.data(), &qp);
    d_ocp_qp_set_q(N, q_stage.data(), &qp);
    d_ocp_qp_set_idxbx(N, idxbx.data(), &qp);
    d_ocp_qp_set_lbx(N, lbx_stage.data(), &qp);
    d_ocp_qp_set_ubx(N, ubx_stage.data(), &qp);

    std::vector<int> idxbxe_all{0, 1, 2, 3};
    d_ocp_qp_set_idxbxe(0, idxbxe_all.data(), &qp);

    auto apply_x0 = [&](const Eigen::Vector<double, NX>& xv)
    {
        Eigen::Vector<double, NX> x0_copy = xv;
        d_ocp_qp_set_lbx(0, x0_copy.data(), &qp);
        d_ocp_qp_set_ubx(0, x0_copy.data(), &qp);
    };
    apply_x0(x0);

    // --- HPIPM solver workspace ----------------------------------------------
    aligned_buffer arg_buf{static_cast<std::size_t>(d_ocp_qp_ipm_arg_memsize(&dim))};
    d_ocp_qp_ipm_arg arg{};
    d_ocp_qp_ipm_arg_create(&dim, &arg, arg_buf.p);
    d_ocp_qp_ipm_arg_set_default(BALANCE, &arg);
    int iter_max = 30;
    d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
    double tol = 1e-6;
    d_ocp_qp_ipm_arg_set_tol_stat(&tol, &arg);
    d_ocp_qp_ipm_arg_set_tol_eq(&tol, &arg);
    d_ocp_qp_ipm_arg_set_tol_ineq(&tol, &arg);
    d_ocp_qp_ipm_arg_set_tol_comp(&tol, &arg);

    aligned_buffer sol_buf{static_cast<std::size_t>(d_ocp_qp_sol_memsize(&dim))};
    d_ocp_qp_sol qp_sol{};
    d_ocp_qp_sol_create(&dim, &qp_sol, sol_buf.p);

    aligned_buffer ws_buf{static_cast<std::size_t>(d_ocp_qp_ipm_ws_memsize(&dim, &arg))};
    d_ocp_qp_ipm_ws ws{};
    d_ocp_qp_ipm_ws_create(&dim, &arg, &ws, ws_buf.p);

    // Warm up
    d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &ws);

    // --- Benchmark ------------------------------------------------------------
    Eigen::Vector<double, NU> hpipm_u0;
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
        .run("hpipm::d_ocp_qp_ipm_solve",
             [&]
             {
                 apply_x0(x0);
                 d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &ws);
                 d_ocp_qp_sol_get_u(0, &qp_sol, hpipm_u0.data());
                 ankerl::nanobench::doNotOptimizeAway(hpipm_u0);
             });

    std::ofstream csv("bench_qp_vs_hpipm.csv");
    bench.render(comma_csv_tpl, csv);
}
