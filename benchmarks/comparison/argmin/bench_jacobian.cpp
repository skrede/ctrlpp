// Finite-difference against analytic constraint partials, on the same
// multiple-shooting problem solved by the same sequential quadratic program.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "arm_accuracy.h"
#include "bench_metrics.h"
#include "shooting_jacobian.h"

#include "nmpc/pendulum.h"
#include "nmpc/double_integrator.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <Eigen/Dense>

#include <span>
#include <memory>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>

namespace
{

namespace arms = ctrlpp::bench::argmin_arms;
namespace problems = ctrlpp::bench::problems::nmpc;

using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

constexpr char const* fd_label = "fd_jacobian";
constexpr char const* analytic_label = "analytic_jacobian";

auto as_span(const Eigen::VectorXd& v) -> std::span<const double>
{
    return {v.data(), static_cast<std::size_t>(v.size())};
}

class slsqp_arm
{
public:
    explicit slsqp_arm(const ctrlpp::nlp_problem<double>& problem)
        : solver{ctrlpp::argmin_settings<double>{}}, result{},
          update{.x0 = Eigen::VectorXd::Zero(problem.n_vars)}
    {
        if(!solver.setup(problem).has_value())
        {
            std::fprintf(stderr, "ctrlpp::argmin_solver setup failed; a measurement taken against a solver that was never set up is meaningless\n");
            std::exit(EXIT_FAILURE);
        }
        solve();
    }

    void solve() { result = solver.solve(update); }

    const ctrlpp::nlp_result<double>& answer() const { return result; }

private:
    ArgminSlsqp solver;
    ctrlpp::nlp_result<double> result;
    ctrlpp::nlp_update<double> update;
};

auto scored(const ctrlpp::nlp_problem<double>& problem, const slsqp_arm& arm, char const* label)
    -> arms::arm_answer
{
    return {max_constraint_violation(problem, as_span(arm.answer().x)), label, arm.answer().x};
}

template <std::size_t NX, std::size_t NU, typename Dynamics>
auto build_problems(const Dynamics& dynamics, int horizon)
    -> std::pair<ctrlpp::nlp_problem<double>, ctrlpp::nlp_problem<double>>
{
    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, NX, NU>>();
    state->x_ref.resize(static_cast<std::size_t>(horizon + 1), ctrlpp::Vector<double, NX>::Zero());
    state->x0 = problems::unit_first_axis_x0<NX>();

    auto differenced = ctrlpp::detail::build_nmpc_problem<double, NX, NU>(
        dynamics, problems::make_nmpc_quadratic_config<NX, NU>(horizon), state);
    auto analytic = differenced;
    analytic.constraint_jacobian = ctrlpp::bench::build_shooting_jacobian<NX, NU>(dynamics, horizon);
    return {differenced, analytic};
}

// nanobench clears its accumulated results whenever the title changes, so a
// per-cell title would leave only the last cell in the rendered file. The cell
// rides in the row name instead and the title names the whole benchmark.
auto cell_tag(const std::string& system_name, std::size_t nx, int horizon) -> std::string
{
    return " " + system_name + " NX=" + std::to_string(nx) + " N=" + std::to_string(horizon);
}

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_jacobian_benchmark(const std::string& system_name, const Dynamics& dynamics, int horizon,
                            ankerl::nanobench::Bench& bench, std::ostream& quality_csv)
{
    auto [differenced, analytic] = build_problems<NX, NU>(dynamics, horizon);
    slsqp_arm fd_arm{differenced};
    slsqp_arm analytic_arm{analytic};
    const std::string tag = cell_tag(system_name, NX, horizon);
    const std::string fd_row = fd_label + tag;
    const std::string analytic_row = analytic_label + tag;
    const std::vector<arms::arm_answer> answers{scored(differenced, fd_arm, fd_row.c_str()),
                                                scored(analytic, analytic_arm, analytic_row.c_str())};
    const double deviation = arms::solution_spread(answers);

    auto run_fd = [&] { fd_arm.solve(); ankerl::nanobench::doNotOptimizeAway(fd_arm.answer().x); };
    auto run_analytic = [&] { analytic_arm.solve(); ankerl::nanobench::doNotOptimizeAway(analytic_arm.answer().x); };

    ctrlpp::bench::run_with_accuracy(bench, arms::deviation_metric, fd_row, deviation, run_fd);
    ctrlpp::bench::run_with_accuracy(bench, arms::deviation_metric, analytic_row, deviation, run_analytic);
    ctrlpp::bench::run_own_criterion_pair(bench, arms::violation_metric, fd_row.c_str(), answers[0].violation,
                                          run_fd, analytic_row.c_str(), answers[1].violation, run_analytic);

    write_quality_csv_row(quality_csv, system_name, "argmin", "slsqp", fd_label, static_cast<int>(NX), horizon,
                          compute_quality_metrics(differenced, fd_arm.answer()));
    write_quality_csv_row(quality_csv, system_name, "argmin", "slsqp", analytic_label, static_cast<int>(NX),
                          horizon, compute_quality_metrics(analytic, analytic_arm.answer()));
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("NLP: finite-difference vs analytic constraint partials")
        .warmup(50)
        .minEpochIterations(50)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    std::ofstream timing_csv("bench_jacobian_timing.csv");
    std::ofstream quality_csv("bench_jacobian_quality.csv");
    write_quality_csv_header(quality_csv);

    const problems::differentiable_double_integrator_2 integrator;
    const problems::differentiable_pendulum_2 pendulum;

    run_jacobian_benchmark<2, 1>("double_integrator", integrator, 10, bench, quality_csv);
    run_jacobian_benchmark<2, 1>("double_integrator", integrator, 20, bench, quality_csv);
    run_jacobian_benchmark<2, 1>("pendulum", pendulum, 10, bench, quality_csv);
    run_jacobian_benchmark<2, 1>("pendulum", pendulum, 5, bench, quality_csv);

    bench.render(ctrlpp::bench::csv_tpl, timing_csv);
}
