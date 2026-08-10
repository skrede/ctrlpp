#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_SCHEDULE_ARMS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_SCHEDULE_ARMS_H

#include "arm_accuracy.h"
#include "bench_metrics.h"

#include "bench_csv.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_problem.h"

#include <argmin/schedule.h>
#include <argmin/solver/options.h>
#include <argmin/solver/cobyla_policy.h>
#include <argmin/solver/kraft_slsqp_policy.h>

#include <Eigen/Dense>

#include <span>
#include <chrono>
#include <string>
#include <cstdio>
#include <cstddef>
#include <utility>
#include <cstdlib>

namespace ctrlpp::bench::argmin_arms
{

template <typename Schedule>
using group_type = argmin::basic_solver_group<Schedule, Eigen::Dynamic,
                                              ctrlpp::argmin_constrained_problem<double>,
                                              argmin::kraft_slsqp_policy<>, argmin::cobyla_policy>;

using baseline_solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

constexpr char const* fallback_label = "fallback_slsqp_cobyla";
constexpr char const* baseline_label = "baseline_slsqp";

struct schedule_probe
{
    arm_answer answer;
    quality_metrics quality;
};

inline auto make_solver_options() -> argmin::solver_options<>
{
    argmin::solver_options<> options;
    options.max_iterations = 500;
    options.set_objective_threshold(1e-6);
    options.set_step_threshold(1e-6);
    return options;
}

inline auto as_span(const Eigen::VectorXd& v) -> std::span<const double>
{
    return {v.data(), static_cast<std::size_t>(v.size())};
}

inline auto converged(argmin::solver_status status) -> bool
{
    return status == argmin::solver_status::converged || status == argmin::solver_status::ftol_reached
        || status == argmin::solver_status::xtol_reached;
}

inline auto bridged(const ctrlpp::nlp_problem<double>& problem) -> ctrlpp::argmin_constrained_problem<double>
{
    ctrlpp::argmin_constrained_problem<double> bridge;
    bridge.partition(problem);
    return bridge;
}

inline void setup_or_exit(baseline_solver& solver, const ctrlpp::nlp_problem<double>& problem)
{
    if(!solver.setup(problem).has_value())
    {
        std::fprintf(stderr, "ctrlpp::argmin_solver setup failed; a measurement taken against a solver that was never set up is meaningless\n");
        std::exit(EXIT_FAILURE);
    }
}

auto quality_of(const ctrlpp::nlp_problem<double>& problem, const auto& result, double wall_ms) -> quality_metrics
{
    return {.objective = result.objective_value,
            .max_constraint_violation = max_constraint_violation(problem, as_span(result.x)),
            .gradient_norm = result.gradient_norm,
            .success = converged(result.status),
            .iterations = static_cast<int>(result.iterations),
            .solve_time_ms = wall_ms};
}

template <typename Schedule, typename Drive>
auto probe_group(const ctrlpp::nlp_problem<double>& problem, const Schedule& schedule, char const* label,
                 Drive&& drive) -> schedule_probe
{
    auto bridge = bridged(problem);
    const Eigen::VectorXd x0 = Eigen::VectorXd::Zero(problem.n_vars);
    auto options = make_solver_options();
    group_type<Schedule> group{bridge, x0, options, schedule};

    const auto started = std::chrono::steady_clock::now();
    const auto result = drive(group, options);
    const double wall_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count();
    return {{max_constraint_violation(problem, as_span(result.x)), label, result.x},
            quality_of(problem, result, wall_ms)};
}

inline auto probe_baseline(const ctrlpp::nlp_problem<double>& problem) -> schedule_probe
{
    baseline_solver solver{ctrlpp::argmin_settings<double>{}};
    setup_or_exit(solver, problem);
    ctrlpp::nlp_update<double> update{.x0 = Eigen::VectorXd::Zero(problem.n_vars)};
    const auto result = solver.solve(update);
    return {{max_constraint_violation(problem, as_span(result.x)), baseline_label, result.x},
            compute_quality_metrics(problem, result)};
}

inline auto fallback_chain() -> argmin::fallback_schedule
{
    argmin::fallback_schedule schedule;
    schedule.stall_threshold = 10;
    return schedule;
}

inline auto time_boxed(std::chrono::microseconds slice) -> argmin::time_boxed_schedule
{
    argmin::time_boxed_schedule schedule;
    schedule.time_slice = slice;
    return schedule;
}

inline void emit_baseline_rows(const ctrlpp::nlp_problem<double>& problem, ankerl::nanobench::Bench& bench,
                               const arm_answer& answer, double spread)
{
    baseline_solver solver{ctrlpp::argmin_settings<double>{}};
    setup_or_exit(solver, problem);
    ctrlpp::nlp_update<double> update{.x0 = Eigen::VectorXd::Zero(problem.n_vars)};
    auto drive = [&] { ankerl::nanobench::doNotOptimizeAway(solver.solve(update)); };

    ctrlpp::bench::run_with_accuracy(bench, spread_metric, answer.label, spread, drive);
    ctrlpp::bench::run_own_criterion_row(bench, violation_metric, answer.label.c_str(), answer.violation, drive);
}

template <typename Schedule, typename Drive>
void emit_group_rows(const ctrlpp::nlp_problem<double>& problem, ankerl::nanobench::Bench& bench,
                     const Schedule& schedule, const arm_answer& answer, double spread, Drive&& drive)
{
    auto bridge = bridged(problem);
    const Eigen::VectorXd x0 = Eigen::VectorXd::Zero(problem.n_vars);
    auto options = make_solver_options();
    auto run = [&]
    {
        group_type<Schedule> group{bridge, x0, options, schedule};
        ankerl::nanobench::doNotOptimizeAway(drive(group, options));
    };

    ctrlpp::bench::run_with_accuracy(bench, spread_metric, answer.label, spread, run);
    ctrlpp::bench::run_own_criterion_row(bench, violation_metric, answer.label.c_str(), answer.violation, run);
}

}

#endif
