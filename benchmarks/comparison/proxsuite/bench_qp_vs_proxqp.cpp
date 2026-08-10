// Competitive benchmark: ctrlpp OSQP wrapper vs proxsuite::proxqp::dense

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "qp/qp_accuracy.h"
#include "qp/qp_reference.h"
#include "qp/dense_mpc_shaped.h"

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <proxsuite/proxqp/dense/dense.hpp>

#include <Eigen/Dense>

#include <cstdio>
#include <fstream>

namespace
{

namespace problems = ctrlpp::bench::problems::qp;

constexpr char const* deviation_metric = "max abs deviation of the two arms' primal solutions x";
constexpr char const* kkt_metric = "relative Karush-Kuhn-Tucker residual of this arm's own primal-dual answer";
constexpr char const* ctrlpp_label = "ctrlpp::osqp_solver::solve";
constexpr char const* proxqp_label = "proxsuite::proxqp::dense::QP::solve";
constexpr double tolerance = 1e-3;
constexpr int iteration_budget = 4000;

// A solver that declined the program performs a fraction of the work and its
// answer scores nothing, so a run that included one would report a figure for a
// problem nobody solved.
bool accepted(ctrlpp::solve_status status)
{
    return status == ctrlpp::solve_status::optimal || status == ctrlpp::solve_status::solved_inaccurate;
}

}

int main(int argc, char** argv)
{
    const problems::dense_program program = problems::make_dense_mpc_program();
    if(!problems::poses_active_constraint(program))
    {
        std::fprintf(stderr, "no constraint row binds at the solution; the rows would not measure constraint handling\n");
        return 1;
    }

    ctrlpp::qp_problem<double> ctrlpp_problem{.P = problems::make_dense_mpc_hessian(),
                                              .q = program.q,
                                              .A = problems::make_dense_mpc_constraint_matrix(),
                                              .l = program.l,
                                              .u = program.u};
    ctrlpp::osqp_solver ctrlpp_solver(tolerance, tolerance, iteration_budget, false, true, true);
    if(!ctrlpp_solver.setup(ctrlpp_problem).has_value())
    {
        std::fprintf(stderr, "ctrlpp::osqp_solver setup failed; a timing measured against a solver that was never set up is meaningless\n");
        return 1;
    }
    ctrlpp::qp_update<double> ctrlpp_update{.q = program.q, .l = program.l, .u = program.u};
    const ctrlpp::qp_result<double> ctrlpp_answer = ctrlpp_solver.solve(ctrlpp_update);

    proxsuite::proxqp::dense::QP<double> proxqp(program.P.cols(), 0, program.A.rows());
    proxqp.settings.eps_abs = tolerance;
    proxqp.settings.eps_rel = tolerance;
    proxqp.settings.verbose = false;
    proxqp.settings.max_iter = iteration_budget;
    proxqp.init(program.P, program.q, std::nullopt, std::nullopt, program.A, program.l, program.u);
    proxqp.solve();
    const Eigen::VectorXd proxqp_x = proxqp.results.x;
    const Eigen::VectorXd proxqp_y = proxqp.results.z;

    if(!accepted(ctrlpp_answer.status) || proxqp.results.info.status != proxsuite::proxqp::QPSolverOutput::PROXQP_SOLVED)
    {
        std::fprintf(stderr, "a solver declined the program; the reported figures would be meaningless\n");
        return 1;
    }

    auto solve_ctrlpp = [&] {
        auto r = ctrlpp_solver.solve(ctrlpp_update);
        ankerl::nanobench::doNotOptimizeAway(r);
    };
    auto solve_proxqp = [&] {
        proxqp.solve();
        auto x = proxqp.results.x;
        ankerl::nanobench::doNotOptimizeAway(x);
    };

    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp vs ProxQP")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::report_accuracy(bench, deviation_metric,
                                   (ctrlpp_answer.x - proxqp_x).cwiseAbs().maxCoeff());
    bench.run(ctrlpp_label, solve_ctrlpp).run(proxqp_label, solve_proxqp);
    ctrlpp::bench::run_own_criterion_pair(
        bench, kkt_metric, ctrlpp_label,
        problems::kkt_relative_residual(program, ctrlpp_answer.x, ctrlpp_answer.y), solve_ctrlpp, proxqp_label,
        problems::kkt_relative_residual(program, proxqp_x, proxqp_y), solve_proxqp);

    std::ofstream csv("bench_qp_vs_proxqp.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
