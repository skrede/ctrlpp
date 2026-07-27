// Churning-active-set QP comparison: ctrlpp OSQP backend vs argmin sparse ADMM.
//
// The sibling bench_qp_vs_osqp_closedloop drives a smoothly moving initial
// condition, so the set of saturated input bounds rotates gently and the warm
// start stays a good predictor -- a stable active-set pattern. THIS bench does
// the opposite: a decorrelated, step-to-step-jumping initial condition drives
// the regulator so a DIFFERENT subset of the input bounds saturates each step.
// The active set churns, the warm-start advantage erodes, and the question is
// whether the sparse-ADMM iteration edge and the polish on/off split survive
// that regime rather than the stable-pattern one.
//
// The churn is MEASURED, not assumed: at every step each input variable is
// classified from the OSQP solution as {at -bound, interior, at +bound}, and the
// per-step fraction of variables whose class changed from the previous step is
// reported as churn_frac (0 = frozen active set, 1 = every input flips class each
// step). A result is only a churn result if churn_frac is high.
//
// Both solvers are handed the identical, deterministic (RNG-free) QP sequence, so
// the only difference measured is solver internals: the same OSQP
// operator-splitting algorithm, argmin's header-only C++ against the vendored
// OSQP C library through ctrlpp::osqp_solver. Both reach the same optimum.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/options/sparse_qp_options.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cstdint>
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

constexpr int kSteps = 60;      // closed-loop steps per cell
constexpr double kULim = 0.3;   // tight input bound: many inputs saturate
constexpr double kAmp = 8.0;    // IC amplitude: drives hard against the bounds
constexpr double kActiveTol = 1e-4; // |u - bound| below this counts as saturated

struct mpc_qp
{
    ctrlpp::qp_problem<double> problem;
    Eigen::VectorXd q;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
    int nx;
    int n_x; // (N+1)*nx: offset of the first input variable in the decision vector
    int nu;
    int horizon;
};

// Same MPC QP as the closed-loop sweep: mildly-coupled stable plant, regulation
// to the origin, dynamics equalities + tight input box bounds. Input-box-only so
// the QP is feasible for any initial condition; the churn lives in which input
// bounds bind, not in problem feasibility.
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
        at.emplace_back(i, i, 1.0); // initial condition, l=u set per step

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

    Eigen::VectorXd q = Eigen::VectorXd::Zero(n_dec);
    return mpc_qp{{.P = P, .q = q, .A = A, .l = l, .u = u}, q, l, u, nx, n_x, nu, N};
}

// Deterministic integer-hash initial condition in [-kAmp, kAmp), decorrelated
// step to step and across components. No RNG -> fully reproducible. Unlike a
// smooth cos driver, consecutive states share almost nothing, so the optimal set
// of saturated inputs churns instead of drifting.
auto churn_ic(int step, int nx) -> Eigen::VectorXd
{
    Eigen::VectorXd x(nx);
    for(int i = 0; i < nx; ++i)
    {
        std::uint32_t h = static_cast<std::uint32_t>(step) * 2654435761u
                        + static_cast<std::uint32_t>(i) * 40503u + 0x9e3779b9u;
        h ^= h >> 15;
        h *= 2246822519u;
        h ^= h >> 13;
        h *= 3266489917u;
        h ^= h >> 16;
        const double f = static_cast<double>(h) / 4294967296.0; // [0,1)
        x(i) = kAmp * (2.0 * f - 1.0);
    }
    return x;
}

// Ternary saturation class of each input variable: -1 at lower bound, +1 at
// upper bound, 0 interior. The active set of the input box is exactly the set of
// non-zero classes.
auto input_classes(const Eigen::VectorXd& sol, int n_x, int nu, int N) -> std::vector<int>
{
    std::vector<int> cls(static_cast<std::size_t>(N * nu), 0);
    for(int k = 0; k < N; ++k)
        for(int i = 0; i < nu; ++i)
        {
            const double v = sol(n_x + k * nu + i);
            const std::size_t idx = static_cast<std::size_t>(k * nu + i);
            if(v <= -kULim + kActiveTol)
                cls[idx] = -1;
            else if(v >= kULim - kActiveTol)
                cls[idx] = 1;
        }
    return cls;
}

// Mean fraction of input variables whose saturation class changes step to step.
auto mean_churn(const std::vector<std::vector<int>>& classes) -> double
{
    if(classes.size() < 2)
        return 0.0;
    double acc = 0.0;
    for(std::size_t s = 1; s < classes.size(); ++s)
    {
        int changed = 0;
        for(std::size_t j = 0; j < classes[s].size(); ++j)
            if(classes[s][j] != classes[s - 1][j])
                ++changed;
        acc += static_cast<double>(changed) / static_cast<double>(classes[s].size());
    }
    return acc / static_cast<double>(classes.size() - 1);
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
    std::ofstream csv("bench_qp_vs_osqp_churning.csv");
    csv << R"("cell","solver","polish","dec","con","iter_min","iter_med","iter_max","iter_mean","elapsed_us","instructions","max_sol_dev_vs_osqp","churn_frac")"
        << "\n";

    ankerl::nanobench::Bench bench;
    bench.title("QP churning active set: ctrlpp OSQP vs argmin sparse ADMM")
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

        for(bool polish : {true, false})
        {
            const char* pol = polish ? "on" : "off";

            // ---- OSQP: setup once, replay the churning IC sequence ----
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
            std::vector<std::vector<int>> classes;
            osqp_iters.reserve(kSteps);
            osqp_sol.reserve(kSteps);
            classes.reserve(kSteps);
            for(int s = 0; s < kSteps; ++s)
            {
                auto xi = churn_ic(s, nx);
                up.l.head(nx) = xi;
                up.u.head(nx) = xi;
                auto r = osqp.solve(up);
                osqp_iters.push_back(r.iterations);
                osqp_sol.push_back(r.x);
                classes.push_back(input_classes(r.x, data.n_x, data.nu, data.horizon));
            }
            const double churn = mean_churn(classes);

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
                auto xi = churn_ic(s, nx);
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
                          auto xi = churn_ic(step_o++ % kSteps, nx);
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
                          auto xi = churn_ic(step_a++ % kSteps, nx);
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
                    << us_of(r) << ',' << instr_of(r) << ',' << dev << ',' << churn << "\n";
            };
            row("ctrlpp::osqp", os, ro, 0.0);
            row("argmin::sparse_admm", as, ra, max_sol_dev);
        }
    }
}
