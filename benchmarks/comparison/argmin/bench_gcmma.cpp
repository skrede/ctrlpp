// bench_gcmma.cpp
//
// Benchmarks GCMMA on the two reference plants (double integrator +
// pendulum) at NX in {2, 4, 8} and horizon N in {10, 20, 30}. NLopt does
// not implement GCMMA; the closest same-paper proxy is CCSAQ, which is the
// affine-plus-quadratic-penalty variant of CCSA (Conservative Convex
// Separable Approximations) from Svanberg (2002), "A class of globally
// convergent optimization methods based on conservative convex separable
// approximations", SIAM J. Optim. 12(2):555-573 — the same paper that
// introduces GCMMA. CCSAQ is therefore "NLopt's recommended globally
// convergent CCSA method", not a strict GCMMA equivalent. Readers
// comparing algorithm families should treat the CCSAQ rows as a proxy
// rather than an equivalence.
//
// Two problem suites per cell:
//
//   Suite A (auglag-wrapped multiple-shooting NMPC, equality continuity):
//     argmin side: NOT INSTANTIATED. The composition
//                  argmin_solver<double, argmin_auglag<argmin_gcmma>>
//                  fails to instantiate against the upstream argmin pin
//                  (probed via `g++ -std=c++23 -c` against the released
//                  milestone/v0.3.0 tag): rho_wval_policy (GCMMA) requires
//                  `constrained<Problem>`, but auglag's inner subproblem
//                  satisfies only `differentiable + bound_constrained`
//                  (it is, by design, the unconstrained AL
//                  re-parameterisation for the inner solver). The two
//                  concepts are mutually exclusive in the current upstream
//                  API. Suite A argmin rows are therefore emitted as
//                  algorithm = "auglag_gcmma" with NaN sentinels and
//                  success=0; the upstream resolution would require
//                  argmin's auglag to expose the zero-constraint stub
//                  interface to inequality-only inner solvers (out of
//                  scope here).
//
//     NLopt side:  nlopt_solver<double> with nlopt_algorithm::auglag_ccsaq
//                  (NLOPT_AUGLAG_EQ outer + NLOPT_LD_CCSAQ inner via
//                  set_local_optimizer; CCSAQ proxy for GCMMA per
//                  Svanberg 2002, see banner above). Equality constraints
//                  are absorbed into the outer AUGLAG_EQ penalty;
//                  inequalities pass through to the inner LD_CCSAQ.
//
//   Suite B (single-shooting unconstrained NLP):
//     argmin side: argmin_solver<double, argmin_gcmma, true>. Note the
//                  `Constrained=true` flag is required even though the
//                  single-shooting problem has zero equality and zero
//                  inequality constraints: upstream rho_wval_policy
//                  static-asserts `constrained<Problem>`, so the bridge
//                  type must be argmin_constrained_problem (which exposes
//                  a zero-row constraint accessor) rather than the
//                  simpler argmin_problem (which exposes none).
//                  Behaviourally equivalent on this NLP because
//                  n_constraints=0; the flag flip is a concept-
//                  satisfaction requirement only.
//     NLopt side:  nlopt_solver<double> with nlopt_algorithm::ccsaq (raw
//                  LD_CCSAQ on the box-constrained, equality-free NLP;
//                  no AUGLAG wrap needed). CCSAQ proxy for GCMMA on
//                  this side as well.
//
// Warm-start configurations per cell: cold, primal_only, curvature.
// `curvature` exercises argmin's GCMMA asymptote-history carry; NLopt has
// no equivalent so the `curvature` NLopt rows are expected to equal the
// `primal_only` NLopt rows per cell. Documented feature, not a bench
// defect.
//
// max_eval / max_time semantics: matched to bench_mma.cpp; see that file's
// header banner. GCMMA's inner/outer split may stall at tight max_eval
// budgets on NX=8 N=30 cells, producing max_iterations status with
// elevated max-violation; this is a known GCMMA limitation on
// constrained multi-shoot, not a bench defect.
//
// Output: bench_gcmma_timing.csv (nanobench mustache) and
// bench_gcmma_quality.csv (bench_metrics.h schema).

#include "bench_metrics.h"
#include "bench_single_shooting.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>
#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

// Common bench knobs applied to every solve so non-converging cells cannot
// stall the harness. max_eval is set above the default 500 so AUGLAG_EQ +
// LD_MMA has room beyond the inner-evaluation accounting ambiguity (per
// 42-02 SUMMARY, the outer counts inner evaluations toward its own
// budget); max_time_sec caps wall-clock per solve so non-converging
// upstream MMA paths (which iterate to the limit on box-only problems
// where MMA's asymptote update has no constraint signal) cannot stall the
// nanobench harness.
constexpr int    bench_max_eval       = 2000;
constexpr double bench_max_time_sec   = 0.5;
constexpr double bench_ftol_rel       = 1e-6;
constexpr double bench_xtol_rel       = 1e-6;
constexpr double bench_constraint_tol = 1e-6;

// ---------------------------------------------------------------------------
// Dynamics (identical to bench_slsqp.cpp / bench_lbfgsb.cpp so suite-A NLopt
// rows are apples-to-apples with the rest of the argmin bench tree).
// ---------------------------------------------------------------------------

constexpr double di2_dt = 0.1;
auto double_integrator_2 = [](const Eigen::Vector2d& x,
                              const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    return Eigen::Vector2d{x(0) + di2_dt * x(1), x(1) + di2_dt * u(0)};
};

auto pendulum_2 = [](const Eigen::Vector2d& x,
                     const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    constexpr double dt = 0.05;
    constexpr double g = 9.81;
    constexpr double l = 1.0;
    double theta = x(0);
    double omega = x(1);
    double alpha = -g / l * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + dt * omega, omega + dt * alpha};
};

constexpr double di4_dt = 0.1;
auto double_integrator_4 = [](const Eigen::Vector4d& x,
                              const Eigen::Vector2d& u) -> Eigen::Vector4d
{
    return Eigen::Vector4d{
        x(0) + di4_dt * x(1),
        x(1) + di4_dt * u(0),
        x(2) + di4_dt * x(3),
        x(3) + di4_dt * u(1)};
};

constexpr double di8_dt = 0.1;
using Vec8 = Eigen::Matrix<double, 8, 1>;
using Vec4 = Eigen::Vector4d;

auto double_integrator_8 = [](const Vec8& x, const Vec4& u) -> Vec8
{
    Vec8 xn;
    xn(0) = x(0) + di8_dt * x(1);
    xn(1) = x(1) + di8_dt * u(0);
    xn(2) = x(2) + di8_dt * x(3);
    xn(3) = x(3) + di8_dt * u(1);
    xn(4) = x(4) + di8_dt * x(5);
    xn(5) = x(5) + di8_dt * u(2);
    xn(6) = x(6) + di8_dt * x(7);
    xn(7) = x(7) + di8_dt * u(3);
    return xn;
};

// ---------------------------------------------------------------------------
// NMPC config factory
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU>
auto make_nmpc_config(int horizon) -> ctrlpp::nmpc_config<double, NX, NU>
{
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix<double, NX, NX>::Identity(),
        .R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1,
    };
}

// ---------------------------------------------------------------------------
// Type aliases
//
// Suite B uses the constrained bridge (Constrained=true) even though the
// single-shooting NLP has n_constraints=0: upstream mma_policy::init
// static-asserts `constrained<Problem>`, and only argmin_constrained_problem
// (the bridge selected by Constrained=true) satisfies that concept.
// ---------------------------------------------------------------------------

using NloptSolver = ctrlpp::nlopt_solver<double>;
using NablaRawGcmma = ctrlpp::argmin_solver<double, ctrlpp::argmin_gcmma, true>;

auto warm_start_label(ctrlpp::warm_start_mode ws) -> std::string
{
    switch(ws)
    {
    case ctrlpp::warm_start_mode::cold:        return "cold";
    case ctrlpp::warm_start_mode::primal_only: return "primal_only";
    case ctrlpp::warm_start_mode::curvature:   return "curvature";
    }
    return "unknown";
}

// Suite-A argmin-side rows are emitted as a static placeholder marking the
// upstream argmin auglag<mma> incompatibility (see banner comment).
void write_argmin_incompatible_row(std::ostream& quality_csv,
                                   const std::string& system,
                                   const std::string& algo_label,
                                   const std::string& ws_label,
                                   int nx,
                                   int horizon)
{
    write_quality_csv_row(quality_csv, system, "argmin", algo_label,
                          ws_label, nx, horizon,
                          quality_metrics{
                              .objective                = std::numeric_limits<double>::quiet_NaN(),
                              .max_constraint_violation = std::numeric_limits<double>::quiet_NaN(),
                              .gradient_norm            = std::numeric_limits<double>::quiet_NaN(),
                              .success                  = false,
                              .iterations               = 0,
                              .solve_time_ms            = 0.0,
                          });
}

// ---------------------------------------------------------------------------
// Suite A: auglag-wrapped multiple-shooting NMPC.
// argmin side is unbuildable upstream (see banner) so only NLopt is timed.
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_suite_a(const std::string& system_name,
                 Dynamics dynamics,
                 int horizon,
                 ankerl::nanobench::Bench& bench,
                 std::ostream& quality_csv,
                 ctrlpp::warm_start_mode ws_mode)
{
    auto config = make_nmpc_config<NX, NU>(horizon);
    auto ws_label = warm_start_label(ws_mode);

    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;

    // MMA paths can stall up to max_time_sec per solve; cut nanobench reps
    // accordingly so the full sweep finishes in single-digit minutes.
    int min_iters    = (NX >= 8) ? 2 : ((NX >= 4) ? 5 : 10);
    int warmup_iters = (NX >= 8) ? 1 : ((NX >= 4) ? 2 : 5);

    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(horizon) + " ws=" + ws_label
               + " suiteA";

    bench.warmup(warmup_iters).minEpochIterations(min_iters).title(title);

    // NLopt AUGLAG_EQ + LD_MMA. NLopt has no warm-start so the same
    // settings run for every ws_mode; the curvature column matches
    // primal_only by construction.
    {
        ctrlpp::nlopt_settings<double> nlopt_cfg{};
        nlopt_cfg.algorithm      = ctrlpp::nlopt_algorithm::auglag_ccsaq;
        nlopt_cfg.ftol_rel       = bench_ftol_rel;
        nlopt_cfg.xtol_rel       = bench_xtol_rel;
        nlopt_cfg.max_eval       = bench_max_eval;
        nlopt_cfg.max_time       = bench_max_time_sec;
        nlopt_cfg.constraint_tol = bench_constraint_tol;

        ctrlpp::nmpc_dynamic<double, NX, NU, NloptSolver, Dynamics> nmpc{
            dynamics, config, NloptSolver{nlopt_cfg}};
        bench.run("nlopt_auglag_ccsaq",
                  [&]
                  {
                      auto u = nmpc.solve(x0);
                      ankerl::nanobench::doNotOptimizeAway(u);
                  });

        ctrlpp::nmpc_dynamic<double, NX, NU, NloptSolver, Dynamics> q{
            dynamics, config, NloptSolver{nlopt_cfg}};
        q.solve(x0);
        auto diag = q.diagnostics();
        auto grad = compute_gradient_norm<double, NX, NU>(q);
        write_quality_csv_row(quality_csv, system_name, "nlopt", "auglag_ccsaq",
                              ws_label, static_cast<int>(NX), horizon,
                              quality_metrics{
                                  .objective                = diag.cost,
                                  .max_constraint_violation = diag.max_constraint_violation,
                                  .gradient_norm            = grad,
                                  .success                  = (diag.status == ctrlpp::solve_status::optimal),
                                  .iterations               = diag.iterations,
                                  .solve_time_ms            = diag.solve_time * 1000.0,
                              });
    }

    // argmin auglag<gcmma> is unbuildable upstream; emit a sentinel row.
    write_argmin_incompatible_row(quality_csv, system_name, "auglag_gcmma",
                                  ws_label, static_cast<int>(NX), horizon);
}

// ---------------------------------------------------------------------------
// Suite B: single-shooting unconstrained NLP, raw MMA on both backends.
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_suite_b(const std::string& system_name,
                 Dynamics dynamics,
                 int horizon,
                 ankerl::nanobench::Bench& bench,
                 std::ostream& quality_csv,
                 ctrlpp::warm_start_mode ws_mode)
{
    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;
    Eigen::Matrix<double, NX, NX> Q = Eigen::Matrix<double, NX, NX>::Identity();
    Eigen::Matrix<double, NU, NU> R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1;

    constexpr double u_min = -10.0;
    constexpr double u_max = 10.0;

    auto problem = bench::build_single_shooting_problem<NX, NU>(
        dynamics, x0, horizon, Q, R, u_min, u_max);

    int n_vars = horizon * static_cast<int>(NU);

    // MMA paths can stall up to max_time_sec per solve; cut nanobench reps
    // accordingly so the full sweep finishes in single-digit minutes.
    int min_iters    = (NX >= 8) ? 2 : ((NX >= 4) ? 5 : 10);
    int warmup_iters = (NX >= 8) ? 1 : ((NX >= 4) ? 2 : 5);

    auto ws_label = warm_start_label(ws_mode);
    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(horizon) + " ws=" + ws_label
               + " suiteB";

    bench.warmup(warmup_iters).minEpochIterations(min_iters).title(title);

    ctrlpp::nlp_update<double> update{};
    update.x0 = Eigen::VectorXd::Zero(n_vars);

    // NLopt raw LD_MMA on the box-constrained, equality-free single-shooting
    // problem. NLopt has no warm-start; settings are identical across
    // ws_mode.
    {
        ctrlpp::nlopt_settings<double> nlopt_cfg{};
        nlopt_cfg.algorithm      = ctrlpp::nlopt_algorithm::ccsaq;
        nlopt_cfg.ftol_rel       = bench_ftol_rel;
        nlopt_cfg.xtol_rel       = bench_xtol_rel;
        nlopt_cfg.max_eval       = bench_max_eval;
        nlopt_cfg.max_time       = bench_max_time_sec;
        nlopt_cfg.constraint_tol = bench_constraint_tol;

        NloptSolver solver{nlopt_cfg};

        bench.run("nlopt_ccsaq",
                  [&]
                  {
                      solver.setup(problem);
                      auto r = solver.solve(update);
                      ankerl::nanobench::doNotOptimizeAway(r);
                  });

        solver.setup(problem);
        auto r = solver.solve(update);
        auto qm = compute_quality_metrics(problem, r);
        write_quality_csv_row(quality_csv, system_name, "nlopt", "ccsaq",
                              ws_label, static_cast<int>(NX), horizon, qm);
    }

    // argmin raw MMA via the constrained bridge; warm_start carried through
    // settings.base; asymptote defaults left untouched (paper-canonical
    // values from argmin_mma_settings).
    {
        ctrlpp::argmin_mma_settings<double> nabla_cfg{};
        nabla_cfg.base.warm_start     = ws_mode;
        nabla_cfg.base.ftol_rel       = bench_ftol_rel;
        nabla_cfg.base.xtol_rel       = bench_xtol_rel;
        nabla_cfg.base.max_eval       = bench_max_eval;
        nabla_cfg.base.max_time       = bench_max_time_sec;
        nabla_cfg.base.constraint_tol = bench_constraint_tol;

        NablaRawGcmma solver{nabla_cfg};

        bench.run("argmin_gcmma",
                  [&]
                  {
                      solver.setup(problem);
                      auto r = solver.solve(update);
                      ankerl::nanobench::doNotOptimizeAway(r);
                  });

        solver.setup(problem);
        auto r = solver.solve(update);
        auto qm = compute_quality_metrics(problem, r);
        write_quality_csv_row(quality_csv, system_name, "argmin", "gcmma",
                              ws_label, static_cast<int>(NX), horizon, qm);
    }
}

// ---------------------------------------------------------------------------
// Convergence reliability: 100 random initial conditions per cell.
// Suite A argmin counted as 0/100 by construction (unbuildable).
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_convergence(const std::string& system_name,
                     Dynamics dynamics,
                     int horizon,
                     std::ostream& quality_csv)
{
    constexpr int num_trials = 100;

    auto config = make_nmpc_config<NX, NU>(horizon);

    ctrlpp::nlopt_settings<double> nlopt_auglag_cfg{};
    nlopt_auglag_cfg.algorithm      = ctrlpp::nlopt_algorithm::auglag_ccsaq;
    nlopt_auglag_cfg.ftol_rel       = bench_ftol_rel;
    nlopt_auglag_cfg.xtol_rel       = bench_xtol_rel;
    nlopt_auglag_cfg.max_eval       = bench_max_eval;
    nlopt_auglag_cfg.max_time       = bench_max_time_sec;
    nlopt_auglag_cfg.constraint_tol = bench_constraint_tol;

    ctrlpp::nlopt_settings<double> nlopt_raw_cfg = nlopt_auglag_cfg;
    nlopt_raw_cfg.algorithm = ctrlpp::nlopt_algorithm::ccsaq;

    ctrlpp::argmin_mma_settings<double> nabla_cfg{};
    nabla_cfg.base.warm_start     = ctrlpp::warm_start_mode::cold;
    nabla_cfg.base.ftol_rel       = bench_ftol_rel;
    nabla_cfg.base.xtol_rel       = bench_xtol_rel;
    nabla_cfg.base.max_eval       = bench_max_eval;
    nabla_cfg.base.max_time       = bench_max_time_sec;
    nabla_cfg.base.constraint_tol = bench_constraint_tol;

    Eigen::Matrix<double, NX, NX> Q = Eigen::Matrix<double, NX, NX>::Identity();
    Eigen::Matrix<double, NU, NU> R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-2.0, 2.0);

    int suite_a_nlopt_ok  = 0;
    int suite_b_argmin_ok = 0;
    int suite_b_nlopt_ok  = 0;

    for(int trial = 0; trial < num_trials; ++trial)
    {
        Eigen::Matrix<double, NX, 1> x0;
        for(std::size_t i = 0; i < NX; ++i)
            x0(static_cast<Eigen::Index>(i)) = dist(rng);

        // Suite A NLopt only (argmin auglag<mma> unbuildable upstream).
        {
            ctrlpp::nmpc_dynamic<double, NX, NU, NloptSolver, Dynamics> nmpc{
                dynamics, config, NloptSolver{nlopt_auglag_cfg}};
            if(nmpc.solve(x0).has_value())
                ++suite_a_nlopt_ok;
        }

        // Suite B
        auto problem = bench::build_single_shooting_problem<NX, NU>(
            dynamics, x0, horizon, Q, R, -10.0, 10.0);
        ctrlpp::nlp_update<double> update{};
        update.x0 = Eigen::VectorXd::Zero(horizon * static_cast<int>(NU));
        {
            NloptSolver solver{nlopt_raw_cfg};
            solver.setup(problem);
            auto r = solver.solve(update);
            if(r.status == ctrlpp::solve_status::optimal)
                ++suite_b_nlopt_ok;
        }
        {
            NablaRawGcmma solver{nabla_cfg};
            solver.setup(problem);
            auto r = solver.solve(update);
            if(r.status == ctrlpp::solve_status::optimal)
                ++suite_b_argmin_ok;
        }
    }

    auto write_rate = [&](char const* solver, char const* algo, int successes)
    {
        double rate = static_cast<double>(successes) / num_trials;
        write_quality_csv_row(quality_csv, system_name, solver, algo,
                              "convergence", static_cast<int>(NX), horizon,
                              quality_metrics{
                                  .objective                = rate,
                                  .max_constraint_violation = 0.0,
                                  .gradient_norm            = 0.0,
                                  .success                  = true,
                                  .iterations               = num_trials,
                                  .solve_time_ms            = 0.0,
                              });
    };

    // argmin auglag<mma> not buildable upstream; record 0/100 explicitly.
    write_rate("argmin", "auglag_gcmma", 0);
    write_rate("nlopt",  "auglag_ccsaq", suite_a_nlopt_ok);
    write_rate("argmin", "gcmma",        suite_b_argmin_ok);
    write_rate("nlopt",  "ccsaq",        suite_b_nlopt_ok);
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_gcmma_timing.csv");
    std::ofstream quality_csv("bench_gcmma_quality.csv");
    write_quality_csv_header(quality_csv);

    // Warm-start trio on representative cells (double_integrator NX=4 N=10
    // and pendulum NX=2 N=10) for both suites; mirrors bench_slsqp.cpp.
    for(auto ws : {ctrlpp::warm_start_mode::cold,
                   ctrlpp::warm_start_mode::primal_only,
                   ctrlpp::warm_start_mode::curvature})
    {
        run_suite_a<4, 2>("double_integrator", double_integrator_4, 10, bench, quality_csv, ws);
        run_suite_a<2, 1>("pendulum",          pendulum_2,          10, bench, quality_csv, ws);
        run_suite_b<4, 2>("double_integrator", double_integrator_4, 10, bench, quality_csv, ws);
        run_suite_b<2, 1>("pendulum",          pendulum_2,          10, bench, quality_csv, ws);
    }

    // Size sweep on double integrator (cold start only, mirrors
    // bench_slsqp.cpp).
    for(int h : {10, 20, 30})
    {
        run_suite_a<2, 1>("double_integrator", double_integrator_2, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
        run_suite_a<4, 2>("double_integrator", double_integrator_4, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
        // NX=8 N=30: GCMMA inner-loop may stall under tight max_eval
        // budgets; report as max_iterations status rather than a bench
        // defect.
        run_suite_a<8, 4>("double_integrator", double_integrator_8, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
        run_suite_b<2, 1>("double_integrator", double_integrator_2, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
        run_suite_b<4, 2>("double_integrator", double_integrator_4, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
        run_suite_b<8, 4>("double_integrator", double_integrator_8, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
    }

    // Pendulum size sweep (NX=2 only; pendulum_4 / pendulum_8 are not in
    // the shared corpus, matching the bench_slsqp.cpp scope).
    for(int h : {20, 30})
    {
        run_suite_a<2, 1>("pendulum", pendulum_2, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
        run_suite_b<2, 1>("pendulum", pendulum_2, h, bench, quality_csv,
                          ctrlpp::warm_start_mode::cold);
    }

    // Convergence reliability (100 random ICs per cell). One representative
    // cell per plant to bound runtime.
    run_convergence<4, 2>("double_integrator", double_integrator_4, 10, quality_csv);
    run_convergence<2, 1>("pendulum",          pendulum_2,          10, quality_csv);

    bench.render(comma_csv_tpl, timing_csv);
    return 0;
}
