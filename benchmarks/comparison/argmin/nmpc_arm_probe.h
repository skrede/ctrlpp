#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_NMPC_ARM_PROBE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_NMPC_ARM_PROBE_H

#include "arm_accuracy.h"
#include "bench_metrics.h"
#include "bench_construct.h"

#include "bench_csv.h"

#include "ctrlpp/nmpc.h"

#include <Eigen/Core>

#include <string>
#include <cstddef>
#include <utility>

namespace ctrlpp::bench::argmin_arms
{

struct nmpc_probe
{
    arm_answer answer;
    quality_metrics quality;
};

/// One untimed solve, from which both the arm's published figures and its
/// quality record are read. The controller is discarded afterwards, so the
/// timed rows below run against a controller in the same fresh state.
template <std::size_t NX, std::size_t NU, typename Solver, typename Dynamics>
auto probe_nmpc_arm(const Dynamics& dynamics, const ctrlpp::nmpc_config<double, NX, NU>& config,
                    const Eigen::Matrix<double, NX, 1>& x0, Solver solver, const std::string& label) -> nmpc_probe
{
    auto controller = ctrlpp::bench::built_or_exit(
        ctrlpp::nmpc_dynamic<double, NX, NU, Solver, Dynamics>::create(dynamics, config, std::move(solver)),
        "arm");
    controller.solve(x0);
    const auto diagnostics = controller.diagnostics();
    return {{compute_constraint_violation(controller), label, controller.last_solution()},
            quality_metrics{.objective = diagnostics.cost,
                            .max_constraint_violation = diagnostics.max_constraint_violation,
                            .gradient_norm = compute_gradient_norm<double, NX, NU>(controller),
                            .success = (diagnostics.status == ctrlpp::solve_status::optimal),
                            .iterations = diagnostics.iterations,
                            .solve_time_ms = diagnostics.solve_time * 1000.0}};
}

template <typename Op>
void emit_variant_rows(ankerl::nanobench::Bench& bench, const arm_answer& answer, double spread, Op&& op)
{
    ctrlpp::bench::run_with_accuracy(bench, spread_metric, answer.label, spread, op);
    ctrlpp::bench::run_own_criterion_row(bench, violation_metric, answer.label.c_str(), answer.violation,
                                         std::forward<Op>(op));
}

}

#endif
