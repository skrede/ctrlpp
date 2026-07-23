// Competitive benchmark: ctrlpp OSQP backend vs argmin sparse ADMM QP solver.
//
// Problem: convex QP with 10 variables, 5 constraints (a representative
// MPC-sized QP), identical data fed to both solvers. Both implement the same
// OSQP operator-splitting algorithm on the canonical form
//   min 1/2 x^T P x + q^T x  s.t.  l <= A x <= u,
// so this measures argmin's header-only C++ implementation against the vendored
// OSQP C library reached through ctrlpp::osqp_solver. Each side is driven
// setup-once / resolve-in-the-hot-loop, the way linear MPC uses it.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/options/sparse_qp_options.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fstream>
#include <vector>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

constexpr int n_vars = 10;
constexpr int n_cons = 5;

auto make_hessian() -> Eigen::SparseMatrix<double>
{
    // Positive-definite diagonal-dominant Hessian.
    Eigen::SparseMatrix<double> H(n_vars, n_vars);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < n_vars; ++i)
    {
        triplets.emplace_back(i, i, 2.0 + 0.1 * i);
        if(i + 1 < n_vars)
        {
            triplets.emplace_back(i, i + 1, 0.1);
            triplets.emplace_back(i + 1, i, 0.1);
        }
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

auto make_constraint_matrix() -> Eigen::SparseMatrix<double>
{
    // Each constraint involves 2 variables.
    Eigen::SparseMatrix<double> A(n_cons, n_vars);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < n_cons; ++i)
    {
        triplets.emplace_back(i, 2 * i, 1.0);
        triplets.emplace_back(i, 2 * i + 1, 0.5);
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

} // namespace

int main()
{
    auto H = make_hessian();
    auto A_con = make_constraint_matrix();
    Eigen::VectorXd q = Eigen::VectorXd::LinSpaced(n_vars, -1.0, 1.0);
    Eigen::VectorXd l = Eigen::VectorXd::Constant(n_cons, -2.0);
    Eigen::VectorXd u = Eigen::VectorXd::Constant(n_cons, 2.0);

    // ---- ctrlpp OSQP backend: setup once, resolve in the loop ----
    ctrlpp::qp_problem<double> problem{
        .P = H,
        .q = q,
        .A = A_con,
        .l = l,
        .u = u};

    ctrlpp::osqp_solver osqp(1e-3, 1e-3, 4000, false, true, true);
    osqp.setup(problem);

    ctrlpp::qp_update<double> osqp_update;
    osqp_update.q = q;
    osqp_update.l = l;
    osqp_update.u = u;
    osqp.solve(osqp_update); // warm up

    // ---- argmin sparse ADMM QP: pose once, resolve in the loop ----
    argmin::sparse_qp_options argmin_opts;
    argmin_opts.eps_abs = 1e-3;
    argmin_opts.eps_rel = 1e-3;
    argmin_opts.max_iterations = 4000;
    argmin_opts.warm_start = true;
    argmin_opts.polish = true;

    argmin::sparse_admm_qp_solver<double> argmin_qp;
    argmin::qp_result<double> argmin_out;
    argmin_qp.solve_into(H, q, A_con, l, u, argmin_out, argmin_opts); // pose + factorize
    argmin_qp.resolve_into(q, l, u, argmin_out, argmin_opts);         // warm up

    // ---- Benchmark: resolve on the frozen factorization, both sides ----
    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp OSQP vs argmin sparse ADMM")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::osqp_solver::solve",
             [&]
             {
                 auto r = osqp.solve(osqp_update);
                 ankerl::nanobench::doNotOptimizeAway(r);
             })
        .run("argmin::sparse_admm_qp_solver::resolve_into",
             [&]
             {
                 argmin_qp.resolve_into(q, l, u, argmin_out, argmin_opts);
                 ankerl::nanobench::doNotOptimizeAway(argmin_out);
             });

    std::ofstream csv("bench_qp_vs_osqp.csv");
    bench.render(comma_csv_tpl, csv);
}
