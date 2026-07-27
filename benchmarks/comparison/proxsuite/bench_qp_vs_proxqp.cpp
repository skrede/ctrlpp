// Competitive benchmark: ctrlpp OSQP wrapper vs proxsuite::proxqp::dense

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "qp/dense_mpc_shaped.h"

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <proxsuite/proxqp/dense/dense.hpp>

#include <Eigen/Dense>

#include <cstdio>
#include <fstream>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

}

int main()
{
    namespace problems = ctrlpp::bench::problems::qp;

    auto H_sparse = problems::make_dense_mpc_hessian();
    auto A_sparse = problems::make_dense_mpc_constraint_matrix();
    auto q = problems::make_dense_mpc_gradient();
    auto lb = problems::make_dense_mpc_lower_bound();
    auto ub = problems::make_dense_mpc_upper_bound();

    Eigen::MatrixXd H = Eigen::MatrixXd(H_sparse);
    Eigen::MatrixXd A = Eigen::MatrixXd(A_sparse);

    constexpr int n = problems::dense_mpc_n_vars;
    constexpr int n_eq = 0;
    constexpr int n_in = problems::dense_mpc_n_cons;

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

    // ---- ProxQP (dense) ----------------------------------------------------
    proxsuite::proxqp::dense::QP<double> proxqp(n, n_eq, n_in);
    proxqp.settings.eps_abs = 1e-3;
    proxqp.settings.eps_rel = 1e-3;
    proxqp.settings.verbose = false;
    proxqp.settings.max_iter = 4000;
    proxqp.init(H, q, std::nullopt, std::nullopt, A, lb, ub);
    proxqp.solve();

    // ---- Benchmark ---------------------------------------------------------
    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp vs ProxQP")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::osqp_solver::solve",
             [&]
             {
                 auto r = ctrlpp_solver.solve(ctrlpp_update);
                 ankerl::nanobench::doNotOptimizeAway(r);
             })
        .run("proxsuite::proxqp::dense::QP::solve",
             [&]
             {
                 proxqp.solve();
                 auto x = proxqp.results.x;
                 ankerl::nanobench::doNotOptimizeAway(x);
             });

    std::ofstream csv("bench_qp_vs_proxqp.csv");
    bench.render(comma_csv_tpl, csv);
}
