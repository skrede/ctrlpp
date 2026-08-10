// Solver schedules on one nonlinear predictive-control problem: a fallback
// chain, the plain sequential quadratic program it falls back from, and three
// time-boxed budgets.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "bench_metrics.h"
#include "schedule_arms.h"

#include "nmpc/double_integrator.h"

#include "ctrlpp/types.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <Eigen/Dense>

#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <fstream>

namespace
{

namespace arms = ctrlpp::bench::argmin_arms;
namespace problems = ctrlpp::bench::problems::nmpc;

struct budget
{
    char const* label;
    std::chrono::microseconds slice;
};

auto budgets() -> std::vector<budget>
{
    return {{"time_boxed_100us", std::chrono::microseconds{100}},
            {"time_boxed_1ms", std::chrono::microseconds{1000}},
            {"time_boxed_10ms", std::chrono::microseconds{10000}}};
}

auto step_to_budget()
{
    return [](auto& group, auto& options) { return group.step_n(500, options); };
}

auto build_problem() -> ctrlpp::nlp_problem<double>
{
    constexpr int horizon = 10;
    auto config = problems::make_nmpc_quadratic_config<4, 2>(horizon);
    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, 4, 2>>();
    state->x_ref.resize(static_cast<std::size_t>(horizon + 1), Eigen::Vector4d::Zero());
    state->x0 = problems::unit_first_axis_x0<4>();
    return ctrlpp::detail::build_nmpc_problem<double, 4, 2>(problems::double_integrator_4, config, state);
}

void write_schedule_row(std::ostream& csv, char const* algorithm, const std::string& label,
                        const quality_metrics& quality)
{
    write_quality_csv_row(csv, "double_integrator", "argmin", algorithm, label, 4, 10, quality);
}

auto probe_budgets(const ctrlpp::nlp_problem<double>& problem) -> std::vector<arms::schedule_probe>
{
    std::vector<arms::schedule_probe> probes;
    for(const budget& budgeted : budgets())
        probes.push_back(
            arms::probe_group(problem, arms::time_boxed(budgeted.slice), budgeted.label, step_to_budget()));
    return probes;
}

void emit_budget_rows(const ctrlpp::nlp_problem<double>& problem, ankerl::nanobench::Bench& bench,
                      std::ostream& quality_csv, const std::vector<arms::schedule_probe>& probes, double spread)
{
    std::size_t index = 0;
    for(const budget& budgeted : budgets())
    {
        arms::emit_group_rows(problem, bench, arms::time_boxed(budgeted.slice), probes[index].answer, spread,
                              step_to_budget());
        write_schedule_row(quality_csv, "schedule", budgeted.label, probes[index].quality);
        ++index;
    }
}

auto answers_of(const arms::schedule_probe& fallback, const arms::schedule_probe& baseline,
                const std::vector<arms::schedule_probe>& boxed) -> std::vector<arms::arm_answer>
{
    std::vector<arms::arm_answer> answers{fallback.answer, baseline.answer};
    for(const arms::schedule_probe& probe : boxed)
        answers.push_back(probe.answer);
    return answers;
}

}

int main(int argc, char** argv)
{
    const ctrlpp::nlp_problem<double> problem = build_problem();
    const auto solve_to_convergence = [](auto& group, auto&) { return group.solve(); };
    const arms::schedule_probe fallback =
        arms::probe_group(problem, arms::fallback_chain(), arms::fallback_label, solve_to_convergence);
    const arms::schedule_probe baseline = arms::probe_baseline(problem);
    const std::vector<arms::schedule_probe> boxed = probe_budgets(problem);
    const double spread = arms::solution_spread(answers_of(fallback, baseline, boxed));

    ankerl::nanobench::Bench bench;
    bench.title("NLP: solver schedules").warmup(20).minEpochIterations(20).performanceCounters(true).relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    std::ofstream timing_csv("bench_schedule_timing.csv");
    std::ofstream quality_csv("bench_schedule_quality.csv");
    write_quality_csv_header(quality_csv);

    arms::emit_group_rows(problem, bench, arms::fallback_chain(), fallback.answer, spread, solve_to_convergence);
    write_schedule_row(quality_csv, "schedule", arms::fallback_label, fallback.quality);
    arms::emit_baseline_rows(problem, bench, baseline.answer, spread);
    write_schedule_row(quality_csv, "slsqp", arms::baseline_label, baseline.quality);
    emit_budget_rows(problem, bench, quality_csv, boxed, spread);

    bench.render(ctrlpp::bench::csv_tpl, timing_csv);
}
