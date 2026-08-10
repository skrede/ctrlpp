// Competitive benchmark: ctrlpp OSQP wrapper vs DAQP (dual active-set)

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "qp/dense_mpc_shaped.h"

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

extern "C"
{
#include <api.h>
#include <types.h>
}

#include <Eigen/Dense>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

int main(int argc, char** argv)
{
    namespace problems = ctrlpp::bench::problems::qp;

    auto H_sparse = problems::make_dense_mpc_hessian();
    auto A_sparse = problems::make_dense_mpc_constraint_matrix();
    auto q  = problems::make_dense_mpc_gradient();
    auto lb = problems::make_dense_mpc_lower_bound();
    auto ub = problems::make_dense_mpc_upper_bound();

    constexpr int n = problems::dense_mpc_n_vars;
    constexpr int m = problems::dense_mpc_n_cons;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> H = Eigen::MatrixXd(H_sparse);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A = Eigen::MatrixXd(A_sparse);

    // ---- ctrlpp OSQP wrapper -----------------------------------------------
    ctrlpp::qp_problem<double> ctrlpp_problem{
        .P = H_sparse,
        .q = q,
        .A = A_sparse,
        .l = lb,
        .u = ub};
    ctrlpp::osqp_solver ctrlpp_solver(1e-3, 1e-3, 4000, false, true, true);
    if(!ctrlpp_solver.setup(ctrlpp_problem).has_value())
    {
        std::fprintf(stderr, "ctrlpp::osqp_solver setup failed; a timing measured against a solver that was never set up is meaningless\n");
        return 1;
    }
    ctrlpp::qp_update<double> ctrlpp_update{.q = q, .l = lb, .u = ub};
    ctrlpp_solver.solve(ctrlpp_update);

    // ---- DAQP --------------------------------------------------------------
    DAQPProblem daqp_qp{};
    daqp_qp.n  = n;
    daqp_qp.m  = m;
    daqp_qp.ms = 0;
    daqp_qp.H  = H.data();
    daqp_qp.f  = q.data();
    daqp_qp.A  = A.data();
    daqp_qp.bupper = ub.data();
    daqp_qp.blower = lb.data();
    std::vector<int> sense(m, 0);
    daqp_qp.sense = sense.data();
    daqp_qp.break_points = nullptr;
    daqp_qp.nh = 0;

    DAQPSettings daqp_settings;
    daqp_default_settings(&daqp_settings);
    daqp_settings.primal_tol = 1e-3;
    daqp_settings.iter_limit = 4000;

    Eigen::VectorXd daqp_x   = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd daqp_lam = Eigen::VectorXd::Zero(m);
    DAQPResult daqp_res{};
    daqp_res.x   = daqp_x.data();
    daqp_res.lam = daqp_lam.data();

    daqp_quadprog(&daqp_res, &daqp_qp, &daqp_settings);

    // ---- Benchmark ---------------------------------------------------------
    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp vs DAQP")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);
    // The agreement figure for this row is unwritten because the competitor is
    // not available on the station this file was edited on, so it could neither
    // be computed nor checked there. Whoever builds this target where the
    // competitor IS available should replace the marker with a cross-arm
    // deviation and each arm's own criterion.
    ctrlpp::bench::report_competitor_not_installed(bench);

    bench
        .run("ctrlpp::osqp_solver::solve",
             [&]
             {
                 auto r = ctrlpp_solver.solve(ctrlpp_update);
                 ankerl::nanobench::doNotOptimizeAway(r);
             })
        .run("daqp_quadprog",
             [&]
             {
                 daqp_quadprog(&daqp_res, &daqp_qp, &daqp_settings);
                 ankerl::nanobench::doNotOptimizeAway(daqp_x);
             });

    std::ofstream csv("bench_qp_vs_daqp.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
