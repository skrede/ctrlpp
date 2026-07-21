#include "bench_metrics.h"

#include "ctrlpp/types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_problem.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <argmin/schedule.h>
#include <argmin/solver/options.h>

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <memory>
#include <string>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

// ---------------------------------------------------------------------------
// Dynamics
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

auto make_solver_options() -> argmin::solver_options<>
{
    argmin::solver_options<> opts;
    opts.max_iterations = 500;
    opts.set_objective_threshold(1e-6);
    opts.set_step_threshold(1e-6);
    return opts;
}

void write_schedule_row(std::ostream& csv,
                        const std::string& schedule_label,
                        double objective,
                        double gradient_norm,
                        double constraint_violation,
                        bool success,
                        int iterations,
                        double solve_time_ms)
{
    quality_metrics qm{
        .objective = objective,
        .max_constraint_violation = constraint_violation,
        .gradient_norm = gradient_norm,
        .success = success,
        .iterations = iterations,
        .solve_time_ms = solve_time_ms,
    };
    write_quality_csv_row(csv, "double_integrator", "argmin", "schedule", schedule_label,
                          4, 10, qm);
}

// ---------------------------------------------------------------------------
// Fallback chains
// ---------------------------------------------------------------------------

void run_fallback_slsqp_cobyla(const ctrlpp::nlp_problem<double>& problem,
                               ankerl::nanobench::Bench& bench,
                               std::ostream& quality_csv)
{
    ctrlpp::argmin_constrained_problem<double> bridge;
    bridge.partition(problem);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(problem.n_vars);
    auto opts = make_solver_options();

    argmin::fallback_schedule sched;
    sched.stall_threshold = 10;

    using group_type = argmin::basic_solver_group<
        argmin::fallback_schedule,
        Eigen::Dynamic,
        ctrlpp::argmin_constrained_problem<double>,
        argmin::kraft_slsqp_policy<>,
        argmin::cobyla_policy>;

    bench.warmup(20).minEpochIterations(20).title("fallback chains")
        .run("fallback_slsqp_cobyla",
             [&]
             {
                 group_type group{bridge, x0, opts, sched};
                 auto result = group.solve();
                 ankerl::nanobench::doNotOptimizeAway(result);
             });

    group_type group{bridge, x0, opts, sched};
    const auto t_start = std::chrono::steady_clock::now();
    auto result = group.solve();
    const auto t_end = std::chrono::steady_clock::now();

    double wall_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    write_schedule_row(quality_csv, "fallback_slsqp_cobyla",
                       result.objective_value, result.gradient_norm,
                       result.constraint_violation,
                       result.status == argmin::solver_status::converged ||
                       result.status == argmin::solver_status::ftol_reached ||
                       result.status == argmin::solver_status::xtol_reached,
                       static_cast<int>(result.iterations), wall_ms);
}

// MMA fallback chains disabled: argmin mma_policy::state_type lacks
// objective_value member required by basic_solver_group. Re-enable when
// argmin fixes mma_policy state_type to include objective_value.
//
// Affected chains: SLSQP->MMA, MMA->COBYLA

// ---------------------------------------------------------------------------
// Baseline
// ---------------------------------------------------------------------------

void run_baseline_slsqp(const ctrlpp::nlp_problem<double>& problem,
                        ankerl::nanobench::Bench& bench,
                        std::ostream& quality_csv)
{
    ctrlpp::argmin_settings<double> cfg{};
    ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp> solver{cfg};
    solver.setup(problem);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(problem.n_vars);

    bench.warmup(20).minEpochIterations(20).title("baseline")
        .run("baseline_slsqp",
             [&]
             {
                 ctrlpp::nlp_update<double> upd{.x0 = x0};
                 auto result = solver.solve(upd);
                 ankerl::nanobench::doNotOptimizeAway(result);
             });

    ctrlpp::nlp_update<double> upd{.x0 = x0};
    auto result = solver.solve(upd);
    auto qm = compute_quality_metrics(problem, result);
    write_quality_csv_row(quality_csv, "double_integrator", "argmin", "slsqp", "baseline_slsqp",
                          4, 10, qm);
}

// ---------------------------------------------------------------------------
// Time-boxed
// ---------------------------------------------------------------------------

void run_time_boxed(const ctrlpp::nlp_problem<double>& problem,
                    std::chrono::microseconds time_slice,
                    const std::string& label,
                    ankerl::nanobench::Bench& bench,
                    std::ostream& quality_csv)
{
    ctrlpp::argmin_constrained_problem<double> bridge;
    bridge.partition(problem);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(problem.n_vars);
    auto opts = make_solver_options();

    argmin::time_boxed_schedule sched;
    sched.time_slice = time_slice;

    using group_type = argmin::basic_solver_group<
        argmin::time_boxed_schedule,
        Eigen::Dynamic,
        ctrlpp::argmin_constrained_problem<double>,
        argmin::kraft_slsqp_policy<>,
        argmin::cobyla_policy>;

    bench.warmup(20).minEpochIterations(20).title("time-boxed")
        .run(label,
             [&]
             {
                 group_type group{bridge, x0, opts, sched};
                 auto result = group.step_n(500, opts);
                 ankerl::nanobench::doNotOptimizeAway(result);
             });

    group_type group{bridge, x0, opts, sched};
    const auto t_start = std::chrono::steady_clock::now();
    auto result = group.step_n(500, opts);
    const auto t_end = std::chrono::steady_clock::now();

    double wall_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    write_schedule_row(quality_csv, label,
                       result.objective_value, result.gradient_norm,
                       result.constraint_violation,
                       result.status == argmin::solver_status::converged ||
                       result.status == argmin::solver_status::ftol_reached ||
                       result.status == argmin::solver_status::xtol_reached,
                       static_cast<int>(result.iterations), wall_ms);
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_schedule_timing.csv");
    std::ofstream quality_csv("bench_schedule_quality.csv");
    write_quality_csv_header(quality_csv);

    constexpr int horizon = 10;
    auto config = ctrlpp::nmpc_config<double, 4, 2>{
        .horizon = horizon,
        .Q = Eigen::Matrix4d::Identity(),
        .R = Eigen::Matrix2d::Identity() * 0.1,
    };

    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, 4, 2>>();
    state->x_ref.resize(static_cast<std::size_t>(horizon + 1), Eigen::Vector4d::Zero());
    Eigen::Vector4d x0_state = Eigen::Vector4d::Zero();
    x0_state(0) = 1.0;
    state->x0 = x0_state;

    // Full NMPC problem (with equality constraints) for constrained chains
    auto problem_full = ctrlpp::detail::build_nmpc_problem<double, 4, 2>(
        double_integrator_4, config, state);

    // Fallback chains (MMA chains disabled -- see comment above)
    run_fallback_slsqp_cobyla(problem_full, bench, quality_csv);

    // Baseline
    run_baseline_slsqp(problem_full, bench, quality_csv);

    // Time-boxed
    using namespace std::chrono_literals;
    run_time_boxed(problem_full, 100us, "time_boxed_100us", bench, quality_csv);
    run_time_boxed(problem_full, 1000us, "time_boxed_1ms", bench, quality_csv);
    run_time_boxed(problem_full, 10000us, "time_boxed_10ms", bench, quality_csv);

    bench.render(comma_csv_tpl, timing_csv);
}
