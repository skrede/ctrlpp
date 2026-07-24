// Per-iteration cost isolation: ctrlpp OSQP backend vs argmin sparse ADMM.
//
// The closed-loop and churning benches measure total per-solve cost, which mixes
// two things: the number of ADMM iterations (a convergence property, contingent on
// the problem and the solver's default scaling/rho) and the cost of each iteration
// (an implementation property, the same-problem trend that transfers). This bench
// separates them.
//
// Method: on a FIXED QP, force each solver to run exactly N iterations -- set the
// tolerances so tight (1e-14) that the stopping test never fires, and cap
// max_iterations at N -- with polishing OFF (polish is a per-solve step that only
// runs on a converged solve, so forcing non-convergence removes it). Sweep N and
// measure the resolve instruction count. Instructions are linear in N:
//
//     instructions(N) = slope * N + intercept
//         slope     = true cost per ADMM iteration (the transferable trend)
//         intercept = fixed per-resolve overhead (warm-start apply, residual and
//                     termination checks, result translation)
//
// ADMM iteration work is fixed-dimension arithmetic independent of the iterate
// values, so the forced-N instruction count is deterministic; the fit is clean.
// The regression itself is done offline from the CSV.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/options/sparse_qp_options.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fstream>
#include <string>
#include <vector>

namespace
{

struct cell
{
    int nx;
    int nu;
    int horizon;
};

const std::vector<cell> cells = {
    {8, 3, 20},   // dec 228
    {12, 4, 30},  // dec 492
};

// Forced iteration counts, evenly spaced and spanning many rho-update intervals so
// any periodic refactorization is amortized into the slope (that is the real
// per-iteration cost a caller pays).
const std::vector<int> Ns = {20, 40, 60, 80, 120, 160, 240, 320};

constexpr double kULim = 0.3;

struct mpc_qp
{
    ctrlpp::qp_problem<double> problem;
    Eigen::VectorXd q, l, u;
    int nx;
};

auto build_mpc_qp(const cell& c) -> mpc_qp
{
    const int nx = c.nx, nu = c.nu, N = c.horizon;
    const int n_x = (N + 1) * nx;
    const int n_dec = n_x + N * nu;
    const int n_dyn = (N + 1) * nx;
    const int n_con = n_dyn + N * nu;

    Eigen::MatrixXd Ad = 0.9 * Eigen::MatrixXd::Identity(nx, nx);
    for(int i = 0; i + 1 < nx; ++i)
        Ad(i, i + 1) = 0.1;
    Eigen::MatrixXd Bd = Eigen::MatrixXd::Zero(nx, nu);
    for(int i = 0; i < nu; ++i)
        Bd(i, i) = 1.0;

    Eigen::SparseMatrix<double> P(n_dec, n_dec);
    std::vector<Eigen::Triplet<double>> pt;
    for(int k = 0; k <= N; ++k)
    {
        const double w = (k == N) ? 10.0 : 1.0;
        for(int i = 0; i < nx; ++i)
            pt.emplace_back(k * nx + i, k * nx + i, w);
    }
    for(int k = 0; k < N; ++k)
        for(int i = 0; i < nu; ++i)
            pt.emplace_back(n_x + k * nu + i, n_x + k * nu + i, 0.1);
    P.setFromTriplets(pt.begin(), pt.end());
    P.makeCompressed();

    Eigen::SparseMatrix<double> A(n_con, n_dec);
    std::vector<Eigen::Triplet<double>> at;
    Eigen::VectorXd l = Eigen::VectorXd::Zero(n_con);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n_con);

    for(int i = 0; i < nx; ++i)
        at.emplace_back(i, i, 1.0);
    for(int k = 0; k < N; ++k)
    {
        const int row = nx + k * nx;
        for(int i = 0; i < nx; ++i)
        {
            at.emplace_back(row + i, (k + 1) * nx + i, 1.0);
            for(int j = 0; j < nx; ++j)
                at.emplace_back(row + i, k * nx + j, -Ad(i, j));
            for(int j = 0; j < nu; ++j)
                at.emplace_back(row + i, n_x + k * nu + j, -Bd(i, j));
        }
    }
    for(int k = 0; k < N; ++k)
    {
        const int row = n_dyn + k * nu;
        for(int i = 0; i < nu; ++i)
        {
            at.emplace_back(row + i, n_x + k * nu + i, 1.0);
            l(row + i) = -kULim;
            u(row + i) = kULim;
        }
    }
    A.setFromTriplets(at.begin(), at.end());
    A.makeCompressed();

    // Fixed, non-trivial initial condition so the forced iterations do real work.
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n_dec);
    for(int i = 0; i < nx; ++i)
    {
        l(i) = 3.0;
        u(i) = 3.0;
    }
    return mpc_qp{{.P = P, .q = q, .A = A, .l = l, .u = u}, q, l, u, nx};
}

} // namespace

int main()
{
    std::ofstream csv("bench_qp_periter.csv");
    csv << R"("cell","solver","forced_N","iters_actual","instructions")" << "\n";

    ankerl::nanobench::Bench bench;
    bench.title("QP per-iteration cost (forced N, polish off)").warmup(10).minEpochIterations(50).performanceCounters(true);
    auto instr_of = [](const ankerl::nanobench::Result& r)
    { return r.median(ankerl::nanobench::Result::Measure::instructions); };

    for(const auto& c : cells)
    {
        auto data = build_mpc_qp(c);
        const std::string tag = "nx=" + std::to_string(c.nx) + " nu=" + std::to_string(c.nu)
                              + " N=" + std::to_string(c.horizon);

        ctrlpp::qp_update<double> up{.q = data.q, .l = data.l, .u = data.u, .warm_x = {}, .warm_y = {}};

        for(int N : Ns)
        {
            // ---- OSQP: eps ~ 0 so it never terminates early; cap at N; no polish; cold ----
            {
                ctrlpp::osqp_solver osqp(1e-14, 1e-14, N, false, false, false);
                osqp.setup(data.problem);
                const int iters = osqp.solve(up).iterations;
                bench.run("osqp " + tag + " N=" + std::to_string(N),
                          [&] { ankerl::nanobench::doNotOptimizeAway(osqp.solve(up)); });
                csv << '"' << tag << "\",\"ctrlpp::osqp\"," << N << ',' << iters << ','
                    << instr_of(bench.results().back()) << "\n";
            }
            // ---- argmin: same forcing ----
            {
                argmin::sparse_qp_options opts;
                opts.eps_abs = 1e-14;
                opts.eps_rel = 1e-14;
                opts.max_iterations = static_cast<std::uint16_t>(N);
                opts.warm_start = false;
                opts.polish = false;
                argmin::sparse_admm_qp_solver<double> aq;
                argmin::qp_result<double> aout;
                aq.solve_into(data.problem.P, data.q, data.problem.A, data.l, data.u, aout, opts);
                const int iters = aout.iterations;
                bench.run("argmin " + tag + " N=" + std::to_string(N),
                          [&] { aq.resolve_into(data.q, data.l, data.u, aout, opts); ankerl::nanobench::doNotOptimizeAway(aout); });
                csv << '"' << tag << "\",\"argmin::sparse_admm\"," << N << ',' << iters << ','
                    << instr_of(bench.results().back()) << "\n";
            }
        }
    }
}
