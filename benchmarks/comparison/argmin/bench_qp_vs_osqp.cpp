// Competitive benchmark: ctrlpp OSQP backend vs argmin sparse ADMM QP solver,
// swept over predictive-control-representative problem sizes.
//
// Both implement the same operator-splitting algorithm, so across the sweep this
// isolates argmin's header-only C++ implementation against the vendored OSQP C
// library reached through ctrlpp::osqp_solver. Each side is driven setup-once /
// resolve-in-the-hot-loop, the way linear predictive control uses it, and both
// are handed identical data.
//
// NOTE (warm-resolve caveat): q/l/u are held fixed across resolves, so after the
// first solve the retained warm start is already optimal and each timed resolve
// converges in near-minimal iterations. This is a best-case, apples-to-apples
// comparison of resolve overhead; the moving and churning initial-condition
// regimes live in the sibling benchmarks of this directory.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "qp/qp_accuracy.h"
#include "qp/sparse_mpc_sweep.h"

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/options/sparse_qp_options.h>

#include <Eigen/Dense>

#include <cstdio>
#include <fstream>
#include <string>

namespace
{

namespace problems = ctrlpp::bench::problems::qp;

constexpr char const* deviation_metric = "max abs deviation of the two arms' primal solutions";
constexpr char const* kkt_metric = "relative Karush-Kuhn-Tucker residual of this arm's own primal-dual answer";
constexpr double tolerance = 1e-3;
constexpr int iteration_budget = 4000;

auto matched_options() -> argmin::sparse_qp_options
{
    argmin::sparse_qp_options options;
    options.eps_abs = tolerance;
    options.eps_rel = tolerance;
    options.max_iterations = iteration_budget;
    options.warm_start = true;
    options.polish = true;
    return options;
}

auto cell_tag(const problems::mpc_sweep_cell& cell, const problems::mpc_sweep_layout& layout) -> std::string
{
    return "nx=" + std::to_string(cell.nx) + " nu=" + std::to_string(cell.nu) + " N="
         + std::to_string(cell.horizon) + " (dec=" + std::to_string(layout.n_dec) + " con="
         + std::to_string(layout.n_con) + ")";
}

class cell_arms
{
public:
    explicit cell_arms(const problems::mpc_sweep_problem& posed)
        : ready{false}, osqp{tolerance, tolerance, iteration_budget, false, true, true}, admm{}, admm_answer{},
          osqp_answer{}, options{matched_options()}, update{.q = posed.q, .l = posed.l, .u = posed.u},
          program{posed}
    {
        ready = osqp.setup(posed.problem).has_value();
        if(!ready)
            return;
        osqp_answer = osqp.solve(update);
        admm.solve_into(posed.problem.P, posed.q, posed.problem.A, posed.l, posed.u, admm_answer, options);
        solve_admm();
    }

    void solve_osqp() { osqp_answer = osqp.solve(update); }

    void solve_admm() { admm.resolve_into(program.q, program.l, program.u, admm_answer, options); }

    bool usable() const { return ready; }

    const Eigen::VectorXd& osqp_primal() const { return osqp_answer.x; }

    const Eigen::VectorXd& osqp_dual() const { return osqp_answer.y; }

    const Eigen::VectorXd& admm_primal() const { return admm_answer.x; }

    const Eigen::VectorXd& admm_dual() const { return admm_answer.y; }

private:
    bool ready;
    ctrlpp::osqp_solver osqp;
    argmin::sparse_admm_qp_solver<double> admm;
    argmin::qp_result<double> admm_answer;
    ctrlpp::qp_result<double> osqp_answer;
    argmin::sparse_qp_options options;
    ctrlpp::qp_update<double> update;
    const problems::mpc_sweep_problem& program;
};

void emit_cell_rows(ankerl::nanobench::Bench& bench, const problems::dense_program& dense, cell_arms& arms,
                    const std::string& tag)
{
    const std::string osqp_label = "ctrlpp::osqp_solver  " + tag;
    const std::string admm_label = "argmin::sparse_admm  " + tag;
    auto solve_osqp = [&] { arms.solve_osqp(); ankerl::nanobench::doNotOptimizeAway(arms.osqp_primal()); };
    auto solve_admm = [&] { arms.solve_admm(); ankerl::nanobench::doNotOptimizeAway(arms.admm_primal()); };

    ctrlpp::bench::report_accuracy(bench, deviation_metric,
                                   (arms.osqp_primal() - arms.admm_primal()).cwiseAbs().maxCoeff());
    bench.run(osqp_label, solve_osqp).run(admm_label, solve_admm);
    ctrlpp::bench::run_own_criterion_pair(
        bench, kkt_metric, osqp_label.c_str(),
        problems::kkt_relative_residual(dense, arms.osqp_primal(), arms.osqp_dual()), solve_osqp,
        admm_label.c_str(), problems::kkt_relative_residual(dense, arms.admm_primal(), arms.admm_dual()),
        solve_admm);
}

bool run_cell(ankerl::nanobench::Bench& bench, const problems::mpc_sweep_cell& cell)
{
    const problems::mpc_sweep_problem posed = problems::build_mpc_sweep_problem(cell);
    cell_arms arms{posed};
    if(!arms.usable())
    {
        std::fprintf(stderr, "ctrlpp::osqp_solver setup failed; a timing measured against a solver that was never set up is meaningless\n");
        return false;
    }
    emit_cell_rows(bench, problems::dense_form(posed), arms, cell_tag(cell, problems::sweep_layout(cell)));
    return true;
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("QP sweep: ctrlpp OSQP vs argmin sparse ADMM (warm resolve)")
        .warmup(20)
        .minEpochIterations(50)
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    for(const problems::mpc_sweep_cell& cell : problems::mpc_sweep_cells())
        if(!run_cell(bench, cell))
            return 1;

    std::ofstream csv("bench_qp_vs_osqp.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
