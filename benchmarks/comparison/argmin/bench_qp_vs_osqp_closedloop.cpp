// Perturbed-initial-condition (closed-loop) QP comparison: ctrlpp OSQP backend
// vs argmin sparse ADMM, over MPC-representative sizes.
//
// The sibling bench_qp_vs_osqp holds q/l/u fixed, so the warm start is already
// optimal and each resolve converges in near-minimal iterations. THIS bench
// instead drives the initial-condition rows of the QP through a moving
// trajectory, so every resolve is a genuinely different problem warm-started
// from the previous solution -- the real per-step MPC workload. argmin named
// the per-step iteration DISTRIBUTION, not just totals, as the single most
// useful signal for its pending per-iteration profiling work, so this captures
// per-resolve iteration counts for both solvers in addition to timing.
//
// Both solvers are handed the identical x_init sequence (identical QPs every
// step), so the only difference measured is solver internals. Same OSQP
// operator-splitting algorithm both sides; argmin's header-only C++ vs the
// vendored OSQP C library via ctrlpp::osqp_solver.
//
// NOTE on argmin's 2026-07-23 polish accept-rule fix (f494b0a, Pareto not
// strict-both): a resolve that lands on an exactly-feasible active-set boundary
// may now polish where it previously declined, which can shift argmin's
// iteration counts here. That is a correctness re-disposition, not a regression
// -- both solvers still reach the same optimum.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/options/sparse_qp_options.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdio>
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
    {2, 1, 10},
    {4, 2, 10},
    {4, 2, 20},
    {8, 3, 20},
    {8, 3, 30},
    {12, 4, 30},
};

constexpr int kSteps = 60; // closed-loop steps per cell

struct mpc_qp
{
    ctrlpp::qp_problem<double> problem;
    Eigen::VectorXd q;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
    int nx;
};

// Same MPC QP as the fixed-data sweep: mildly-coupled stable plant, regulation
// to the origin, n_dec = (N+1)nx + N*nu, dynamics equalities + input bounds.
auto build_mpc_qp(const cell& c) -> mpc_qp
{
    const int nx = c.nx;
    const int nu = c.nu;
    const int N = c.horizon;
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

    const double q_weight = 1.0, r_weight = 0.1, qf_weight = 10.0;

    Eigen::SparseMatrix<double> P(n_dec, n_dec);
    std::vector<Eigen::Triplet<double>> pt;
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

    Eigen::SparseMatrix<double> A(n_con, n_dec);
    std::vector<Eigen::Triplet<double>> at;
    Eigen::VectorXd l = Eigen::VectorXd::Zero(n_con);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n_con);

    for(int i = 0; i < nx; ++i)
    {
        at.emplace_back(i, i, 1.0); // initial condition, l=u set per step
    }
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
    // Tight input bounds so the regulator saturates against them: the QP has a
    // moving active set step to step, which is where warm-resolve iteration
    // counts vary and where argmin's polish accept-rule (f494b0a) actually
    // engages. A loose bound leaves an easy interior QP that converges inside a
    // single termination-check interval every step (no distribution to see).
    const double u_lim = 0.5;
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

    Eigen::VectorXd q = Eigen::VectorXd::Zero(n_dec);
    return mpc_qp{{.P = P, .q = q, .A = A, .l = l, .u = u}, q, l, u, nx};
}

// Deterministic moving initial condition: bounded, changes every step, so each
// resolve is a distinct QP. No RNG -> fully reproducible.
auto x_init_at(int step, int nx) -> Eigen::VectorXd
{
    // Large-amplitude, fast-varying initial state so the regulator drives hard
    // against the tight input bounds and the warm start from the previous step
    // is genuinely displaced -- both push resolves past a single termination
    // interval on the harder steps.
    Eigen::VectorXd x(nx);
    for(int i = 0; i < nx; ++i)
        x(i) = 5.0 * std::cos(0.8 * step + 0.7 * i);
    return x;
}

struct iter_stats
{
    int lo, med, hi;
    double mean;
};

auto summarize(std::vector<int> v) -> iter_stats
{
    std::sort(v.begin(), v.end());
    double sum = 0;
    for(int x : v)
        sum += x;
    return {v.front(), v[v.size() / 2], v.back(), sum / static_cast<double>(v.size())};
}

} // namespace

int main()
{
    std::ofstream csv("bench_qp_vs_osqp_closedloop.csv");
    csv << R"("cell","solver","polish","dec","con","iter_min","iter_med","iter_max","iter_mean","elapsed_us","instructions","max_sol_dev_vs_osqp")"
        << "\n";

    ankerl::nanobench::Bench bench;
    bench.title("QP closed-loop (perturbed IC): ctrlpp OSQP vs argmin sparse ADMM")
        .warmup(10)
        .minEpochIterations(static_cast<uint64_t>(kSteps))
        .performanceCounters(true);

    auto us_of = [](const ankerl::nanobench::Result& r)
    { return r.median(ankerl::nanobench::Result::Measure::elapsed) * 1e6; };
    auto instr_of = [](const ankerl::nanobench::Result& r)
    { return r.median(ankerl::nanobench::Result::Measure::instructions); };

    for(const auto& c : cells)
    {
        auto data = build_mpc_qp(c);
        const int nx = data.nx;
        const int n_dec = static_cast<int>(data.q.size());
        const int n_con = static_cast<int>(data.l.size());
        const std::string tag = "nx=" + std::to_string(c.nx) + " nu=" + std::to_string(c.nu)
                              + " N=" + std::to_string(c.horizon);

        // argmin's 2026-07-23 intel: on a moving IC the polish step can dominate
        // its per-step cost (a full symbolic-analysis + factorization of the
        // reduced KKT every resolve). OSQP ships polish off by default; argmin
        // on. We measure BOTH polish states so the ADMM kernel and the polish
        // are separable, as argmin asked.
        for(bool polish : {true, false})
        {
            const char* pol = polish ? "on" : "off";

            // ---- OSQP: setup once (polish fixed at construction), replay ----
            ctrlpp::osqp_solver osqp(1e-3, 1e-3, 4000, false, true, polish);
            if(!osqp.setup(data.problem).has_value())
            {
                std::fprintf(stderr, "ctrlpp::osqp_solver setup failed; a timing measured against a solver that was never set up is meaningless\n");
                return 1;
            }
            ctrlpp::qp_update<double> up;
            up.q = data.q;
            up.l = data.l;
            up.u = data.u;

            std::vector<int> osqp_iters;
            std::vector<Eigen::VectorXd> osqp_sol;
            osqp_iters.reserve(kSteps);
            osqp_sol.reserve(kSteps);
            for(int s = 0; s < kSteps; ++s)
            {
                auto xi = x_init_at(s, nx);
                up.l.head(nx) = xi;
                up.u.head(nx) = xi;
                auto r = osqp.solve(up);
                osqp_iters.push_back(r.iterations);
                osqp_sol.push_back(r.x);
            }

            // ---- argmin: pose once, replay (polish flag in opts) ----
            argmin::sparse_qp_options opts;
            opts.eps_abs = 1e-3;
            opts.eps_rel = 1e-3;
            opts.max_iterations = 4000;
            opts.warm_start = true;
            opts.polish = polish;
            argmin::sparse_admm_qp_solver<double> aq;
            argmin::qp_result<double> aout;
            aq.solve_into(data.problem.P, data.q, data.problem.A, data.l, data.u, aout, opts);

            Eigen::VectorXd al = data.l, au = data.u;
            std::vector<int> argmin_iters;
            argmin_iters.reserve(kSteps);
            double max_sol_dev = 0.0;
            for(int s = 0; s < kSteps; ++s)
            {
                auto xi = x_init_at(s, nx);
                al.head(nx) = xi;
                au.head(nx) = xi;
                aq.resolve_into(data.q, al, au, aout, opts);
                argmin_iters.push_back(aout.iterations);
                max_sol_dev = std::max(max_sol_dev, (aout.x - osqp_sol[static_cast<std::size_t>(s)]).cwiseAbs().maxCoeff());
            }

            auto os = summarize(osqp_iters);
            auto as = summarize(argmin_iters);

            int step_o = 0;
            bench.run("osqp  " + tag + " polish=" + pol,
                      [&]
                      {
                          auto xi = x_init_at(step_o++ % kSteps, nx);
                          up.l.head(nx) = xi;
                          up.u.head(nx) = xi;
                          auto r = osqp.solve(up);
                          ankerl::nanobench::doNotOptimizeAway(r);
                      });
            const auto ro = bench.results().back();

            int step_a = 0;
            bench.run("argmin  " + tag + " polish=" + pol,
                      [&]
                      {
                          auto xi = x_init_at(step_a++ % kSteps, nx);
                          al.head(nx) = xi;
                          au.head(nx) = xi;
                          aq.resolve_into(data.q, al, au, aout, opts);
                          ankerl::nanobench::doNotOptimizeAway(aout);
                      });
            const auto ra = bench.results().back();

            auto row = [&](const char* solver, const iter_stats& st,
                           const ankerl::nanobench::Result& r, double dev)
            {
                csv << '"' << tag << "\",\"" << solver << "\",\"" << pol << "\"," << n_dec << ',' << n_con << ','
                    << st.lo << ',' << st.med << ',' << st.hi << ',' << st.mean << ','
                    << us_of(r) << ',' << instr_of(r) << ',' << dev << "\n";
            };
            row("ctrlpp::osqp", os, ro, 0.0);
            row("argmin::sparse_admm", as, ra, max_sol_dev);
        }
    }
}
