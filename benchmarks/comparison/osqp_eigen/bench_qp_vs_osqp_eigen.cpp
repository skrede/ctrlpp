// Competitive benchmark: ctrlpp OSQP wrapper vs osqp-eigen
// Problem: Convex QP with 10 variables, 5 constraints (representative MPC-sized QP)
// Both ultimately call OSQP v1.0.0 underneath -- this measures wrapper overhead.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <OsqpEigen/OsqpEigen.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstdio>
#include <fstream>

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
    // Positive-definite diagonal-dominant Hessian
    Eigen::SparseMatrix<double> H(n_vars, n_vars);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < n_vars; ++i) {
        triplets.emplace_back(i, i, 2.0 + 0.1 * i);
        if (i + 1 < n_vars) {
            triplets.emplace_back(i, i + 1, 0.1);
            triplets.emplace_back(i + 1, i, 0.1);
        }
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

auto make_constraint_matrix() -> Eigen::SparseMatrix<double>
{
    // Sparse constraint matrix: each constraint involves 2-3 variables
    Eigen::SparseMatrix<double> A(n_cons, n_vars);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < n_cons; ++i) {
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

    // ---- ctrlpp OSQP wrapper setup ----
    ctrlpp::qp_problem<double> problem{
        .P = H,
        .q = q,
        .A = A_con,
        .l = l,
        .u = u};

    ctrlpp::osqp_solver ctrlpp_solver(1e-3, 1e-3, 4000, false, true, true);
    if(!ctrlpp_solver.setup(problem).has_value())
    {
        std::fprintf(stderr, "ctrlpp::osqp_solver setup failed; a timing measured against a solver that was never set up is meaningless\n");
        return 1;
    }

    ctrlpp::qp_update<double> ctrlpp_update;
    ctrlpp_update.q = q;
    ctrlpp_update.l = l;
    ctrlpp_update.u = u;

    // Warm up
    ctrlpp_solver.solve(ctrlpp_update);

    // ---- osqp-eigen setup ----
    OsqpEigen::Solver osqp_eigen_solver;
    osqp_eigen_solver.settings()->setVerbosity(false);
    osqp_eigen_solver.settings()->setWarmStart(true);
    osqp_eigen_solver.settings()->setAbsoluteTolerance(1e-3);
    osqp_eigen_solver.settings()->setRelativeTolerance(1e-3);
    osqp_eigen_solver.settings()->setMaxIteration(4000);
    osqp_eigen_solver.settings()->setPolish(true);
    osqp_eigen_solver.data()->setNumberOfVariables(n_vars);
    osqp_eigen_solver.data()->setNumberOfConstraints(n_cons);
    osqp_eigen_solver.data()->setHessianMatrix(H);
    osqp_eigen_solver.data()->setGradient(q);
    osqp_eigen_solver.data()->setLinearConstraintsMatrix(A_con);
    osqp_eigen_solver.data()->setLowerBound(l);
    osqp_eigen_solver.data()->setUpperBound(u);
    osqp_eigen_solver.initSolver();

    // Warm up
    osqp_eigen_solver.solve();

    // ---- Benchmark ----
    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp vs osqp-eigen")
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
        .run("OsqpEigen::Solver::solve",
             [&]
             {
                 osqp_eigen_solver.solve();
                 auto sol = osqp_eigen_solver.getSolution();
                 ankerl::nanobench::doNotOptimizeAway(sol);
             });

    std::ofstream csv("bench_qp_vs_osqp_eigen.csv");
    bench.render(comma_csv_tpl, csv);
}
