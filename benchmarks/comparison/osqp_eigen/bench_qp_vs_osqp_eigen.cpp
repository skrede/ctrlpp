// Competitive benchmark: ctrlpp OSQP wrapper vs osqp-eigen.
// Both ultimately call OSQP v1.0.0 underneath -- this measures wrapper overhead.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "qp/qp_accuracy.h"
#include "qp/qp_reference.h"
#include "qp/dense_mpc_shaped.h"

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <OsqpEigen/OsqpEigen.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstdio>
#include <fstream>

namespace
{

namespace problems = ctrlpp::bench::problems::qp;

constexpr char const* deviation_metric = "max abs deviation of the two arms' primal solutions x";
constexpr char const* kkt_metric = "relative Karush-Kuhn-Tucker residual of this arm's own primal-dual answer";
constexpr char const* ctrlpp_label = "ctrlpp::osqp_solver::solve";
constexpr char const* osqp_eigen_label = "OsqpEigen::Solver::solve";
constexpr double tolerance = 1e-3;
constexpr int iteration_budget = 4000;

class ctrlpp_arm
{
public:
    explicit ctrlpp_arm(const problems::dense_program& program)
        : ready{false}, solver{tolerance, tolerance, iteration_budget, false, true, true}, answer{},
          update{.q = program.q, .l = program.l, .u = program.u},
          problem{.P = problems::make_dense_mpc_hessian(),
                  .q = program.q,
                  .A = problems::make_dense_mpc_constraint_matrix(),
                  .l = program.l,
                  .u = program.u}
    {
        ready = solver.setup(problem).has_value();
        if(ready)
            answer = solver.solve(update);
    }

    void solve() { answer = solver.solve(update); }

    const Eigen::VectorXd& primal() const { return answer.x; }

    const Eigen::VectorXd& dual() const { return answer.y; }

    bool usable() const
    {
        return ready && (answer.status == ctrlpp::solve_status::optimal
                         || answer.status == ctrlpp::solve_status::solved_inaccurate);
    }

private:
    bool ready;
    ctrlpp::osqp_solver solver;
    ctrlpp::qp_result<double> answer;
    ctrlpp::qp_update<double> update;
    ctrlpp::qp_problem<double> problem;
};

// osqp-eigen's data setters take non-const references, so the operands live
// beside the solver rather than in the caller's frame.
class osqp_eigen_arm
{
public:
    explicit osqp_eigen_arm(const problems::dense_program& program)
        : lower{program.l}, upper{program.u}, gradient{program.q},
          hessian{problems::make_dense_mpc_hessian()},
          constraints{problems::make_dense_mpc_constraint_matrix()}
    {
        configure(program);
        solver.solve();
    }

    void solve() { solver.solve(); }

    const Eigen::VectorXd& primal() { return solver.getSolution(); }

    const Eigen::VectorXd& dual() { return solver.getDualSolution(); }

private:
    Eigen::VectorXd lower;
    Eigen::VectorXd upper;
    OsqpEigen::Solver solver;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> hessian;
    Eigen::SparseMatrix<double> constraints;

    void configure(const problems::dense_program& program)
    {
        solver.settings()->setVerbosity(false);
        solver.settings()->setWarmStart(true);
        solver.settings()->setAbsoluteTolerance(tolerance);
        solver.settings()->setRelativeTolerance(tolerance);
        solver.settings()->setMaxIteration(iteration_budget);
        solver.settings()->setPolish(true);
        solver.data()->setNumberOfVariables(static_cast<int>(program.P.cols()));
        solver.data()->setNumberOfConstraints(static_cast<int>(program.A.rows()));
        solver.data()->setHessianMatrix(hessian);
        solver.data()->setGradient(gradient);
        solver.data()->setLinearConstraintsMatrix(constraints);
        solver.data()->setLowerBound(lower);
        solver.data()->setUpperBound(upper);
        solver.initSolver();
    }
};

void emit_rows(ankerl::nanobench::Bench& bench, const problems::dense_program& program, ctrlpp_arm& mine,
               osqp_eigen_arm& theirs)
{
    auto solve_mine = [&] { mine.solve(); ankerl::nanobench::doNotOptimizeAway(mine.primal()); };
    auto solve_theirs = [&] { theirs.solve(); ankerl::nanobench::doNotOptimizeAway(theirs.primal()); };

    ctrlpp::bench::report_accuracy(bench, deviation_metric,
                                   (mine.primal() - theirs.primal()).cwiseAbs().maxCoeff());
    bench.run(ctrlpp_label, solve_mine).run(osqp_eigen_label, solve_theirs);
    ctrlpp::bench::run_own_criterion_pair(
        bench, kkt_metric, ctrlpp_label, problems::kkt_relative_residual(program, mine.primal(), mine.dual()),
        solve_mine, osqp_eigen_label, problems::kkt_relative_residual(program, theirs.primal(), theirs.dual()),
        solve_theirs);
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

    ctrlpp_arm mine{program};
    if(!mine.usable())
    {
        std::fprintf(stderr, "ctrlpp::osqp_solver declined the program; the reported figures would be meaningless\n");
        return 1;
    }
    osqp_eigen_arm theirs{program};

    ankerl::nanobench::Bench bench;
    bench.title("QP: ctrlpp vs osqp-eigen").warmup(50).minEpochIterations(100).performanceCounters(true).relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);
    emit_rows(bench, program, mine, theirs);

    std::ofstream csv("bench_qp_vs_osqp_eigen.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
