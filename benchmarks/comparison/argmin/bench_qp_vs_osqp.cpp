// Competitive benchmark: ctrlpp OSQP backend vs argmin sparse ADMM QP solver,
// swept over MPC-representative problem sizes.
//
// Each cell builds a genuine sparse linear-MPC QP in canonical OSQP form
//   min 1/2 z^T P z + q^T z  s.t.  l <= A z <= u,
// with the decision vector z = [x_0..x_N, u_0..u_{N-1}] dimensioned exactly as
// ctrlpp's mpc<> lays it out: n_dec = (N+1)*nx + N*nu, and constraints
// n_con = (N+1)*nx (initial condition + dynamics equalities) + N*nu (input
// bounds). Regulation of a mildly-coupled stable plant from a nonzero initial
// state to the origin, so every cell is feasible and both solvers converge to
// the same optimum. Identical data is handed to both.
//
// Both implement the same OSQP operator-splitting algorithm, so across the
// sweep this isolates argmin's header-only C++ implementation against the
// vendored OSQP C library reached through ctrlpp::osqp_solver. Each side is
// driven setup-once / resolve-in-the-hot-loop, the way linear MPC uses it.
//
// NOTE (warm-resolve caveat): q/l/u are held fixed across resolves, so after the
// first solve the retained warm start is already optimal and each timed resolve
// converges in near-minimal iterations. This is a best-case, apples-to-apples
// comparison of resolve overhead; a perturbed-initial-condition variant (which
// would time real per-step MPC work) is a deliberate future refinement.

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

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

struct cell
{
    int nx;
    int nu;
    int horizon;
};

// MPC-representative sweep: tiny-embedded through mid-scale, matching the
// NX=2/4/8/12, N=10/20/30 cells the argmin SQP/NMPC benches already use.
const std::vector<cell> cells = {
    {2, 1, 10},
    {4, 2, 10},
    {4, 2, 20},
    {8, 3, 20},
    {8, 3, 30},
    {12, 4, 30},
};

// A mildly-coupled stable plant: A = 0.9 I + 0.1 superdiagonal shift (spectral
// radius < 1), B = leading nu columns of the identity. Regulated to the origin
// from x_init = 1, with generous input bounds so the QP is always feasible.
struct mpc_qp
{
    ctrlpp::qp_problem<double> problem;
    Eigen::VectorXd q;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
};

auto build_mpc_qp(const cell& c) -> mpc_qp
{
    const int nx = c.nx;
    const int nu = c.nu;
    const int N = c.horizon;

    const int n_x = (N + 1) * nx;
    const int n_dec = n_x + N * nu;
    const int n_dyn = (N + 1) * nx;    // x_0 = x_init, then N transitions
    const int n_con = n_dyn + N * nu;  // + input bounds

    // Plant.
    Eigen::MatrixXd Ad = 0.9 * Eigen::MatrixXd::Identity(nx, nx);
    for(int i = 0; i + 1 < nx; ++i)
        Ad(i, i + 1) = 0.1;
    Eigen::MatrixXd Bd = Eigen::MatrixXd::Zero(nx, nu);
    for(int i = 0; i < nu; ++i)
        Bd(i, i) = 1.0;

    const double q_weight = 1.0;
    const double r_weight = 0.1;
    const double qf_weight = 10.0;

    // Cost P = blkdiag(Q on x_0..x_{N-1}, Qf on x_N, R on u_0..u_{N-1}).
    Eigen::SparseMatrix<double> P(n_dec, n_dec);
    std::vector<Eigen::Triplet<double>> pt;
    pt.reserve(static_cast<std::size_t>(n_dec));
    for(int k = 0; k <= N; ++k)
    {
        const double w = (k == N) ? qf_weight : q_weight;
        for(int i = 0; i < nx; ++i)
            pt.emplace_back(k * nx + i, k * nx + i, w);
    }
    for(int k = 0; k < N; ++k)
        for(int i = 0; i < nu; ++i)
            pt.emplace_back(n_x + k * nu + i, n_x + k * nu + i, r_weight);
    P.setFromTriplets(pt.begin(), pt.end());
    P.makeCompressed();

    // Constraints.
    Eigen::SparseMatrix<double> A(n_con, n_dec);
    std::vector<Eigen::Triplet<double>> at;
    Eigen::VectorXd l = Eigen::VectorXd::Zero(n_con);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n_con);

    const double x_init = 1.0;
    // Initial condition: x_0 = x_init (rows 0..nx).
    for(int i = 0; i < nx; ++i)
    {
        at.emplace_back(i, i, 1.0);
        l(i) = x_init;
        u(i) = x_init;
    }
    // Dynamics: x_{k+1} - Ad x_k - Bd u_k = 0 (rows nx + k*nx).
    for(int k = 0; k < N; ++k)
    {
        const int row = nx + k * nx;
        for(int i = 0; i < nx; ++i)
        {
            at.emplace_back(row + i, (k + 1) * nx + i, 1.0);            // x_{k+1}
            for(int j = 0; j < nx; ++j)
                at.emplace_back(row + i, k * nx + j, -Ad(i, j));        // -Ad x_k
            for(int j = 0; j < nu; ++j)
                at.emplace_back(row + i, n_x + k * nu + j, -Bd(i, j));  // -Bd u_k
        }
        // l = u = 0 already.
    }
    // Input bounds: -u_lim <= u_k <= u_lim (rows n_dyn + k*nu).
    const double u_lim = 10.0;
    for(int k = 0; k < N; ++k)
    {
        const int row = n_dyn + k * nu;
        for(int i = 0; i < nu; ++i)
        {
            at.emplace_back(row + i, n_x + k * nu + i, 1.0);
            l(row + i) = -u_lim;
            u(row + i) = u_lim;
        }
    }
    A.setFromTriplets(at.begin(), at.end());
    A.makeCompressed();

    Eigen::VectorXd q = Eigen::VectorXd::Zero(n_dec); // regulation to the origin

    return mpc_qp{
        .problem = {.P = P, .q = q, .A = A, .l = l, .u = u},
        .q = q,
        .l = l,
        .u = u};
}

} // namespace

int main()
{
    ankerl::nanobench::Bench bench;
    bench.title("QP sweep: ctrlpp OSQP vs argmin sparse ADMM (warm resolve)")
        .warmup(20)
        .minEpochIterations(50)
        .performanceCounters(true);

    for(const auto& c : cells)
    {
        auto data = build_mpc_qp(c);
        const int n_dec = static_cast<int>(data.q.size());
        const int n_con = static_cast<int>(data.l.size());
        const std::string tag = "nx=" + std::to_string(c.nx) + " nu=" + std::to_string(c.nu)
                              + " N=" + std::to_string(c.horizon)
                              + " (dec=" + std::to_string(n_dec) + " con=" + std::to_string(n_con) + ")";

        // ---- ctrlpp OSQP backend: setup once, resolve in the loop ----
        ctrlpp::osqp_solver osqp(1e-3, 1e-3, 4000, false, true, true);
        osqp.setup(data.problem);
        ctrlpp::qp_update<double> osqp_update;
        osqp_update.q = data.q;
        osqp_update.l = data.l;
        osqp_update.u = data.u;
        osqp.solve(osqp_update); // warm up

        // ---- argmin sparse ADMM QP: pose once, resolve in the loop ----
        argmin::sparse_qp_options opts;
        opts.eps_abs = 1e-3;
        opts.eps_rel = 1e-3;
        opts.max_iterations = 4000;
        opts.warm_start = true;
        opts.polish = true;
        argmin::sparse_admm_qp_solver<double> argmin_qp;
        argmin::qp_result<double> argmin_out;
        argmin_qp.solve_into(data.problem.P, data.q, data.problem.A, data.l, data.u, argmin_out, opts);
        argmin_qp.resolve_into(data.q, data.l, data.u, argmin_out, opts); // warm up

        bench.run("ctrlpp::osqp_solver  " + tag,
                  [&]
                  {
                      auto r = osqp.solve(osqp_update);
                      ankerl::nanobench::doNotOptimizeAway(r);
                  });
        bench.run("argmin::sparse_admm  " + tag,
                  [&]
                  {
                      argmin_qp.resolve_into(data.q, data.l, data.u, argmin_out, opts);
                      ankerl::nanobench::doNotOptimizeAway(argmin_out);
                  });
    }

    std::ofstream csv("bench_qp_vs_osqp.csv");
    bench.render(comma_csv_tpl, csv);
}
