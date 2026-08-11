// The sequential-quadratic-programming variants against each other and against
// the reference nonlinear-programming library, on the same predictive-control
// problem across problem sizes.
//
// This target and bench_slsqp are the two halves of one comparison: they pose
// the identical problems from the identical initial states over the same variant
// list, and differ only in warm-start mode. The argmin arms here run at the
// settings default, which retains curvature, so this is the warm half and
// bench_slsqp's size sweep is the cold one. The reference library carries no
// such state and so is timed from scratch in both.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"

#include "arm_accuracy.h"
#include "bench_metrics.h"
#include "bench_construct.h"
#include "nmpc_arm_probe.h"
#include "convergence_rate.h"

#include "nmpc/pendulum.h"
#include "nmpc/double_integrator.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <cstddef>
#include <fstream>

namespace
{

namespace arms = ctrlpp::bench::argmin_arms;
namespace problems = ctrlpp::bench::problems::nmpc;

using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;
using ArgminNwSqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp>;
using ArgminFilterSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_slsqp>;
using ArgminFilterNwSqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_filter_nw_sqp>;
using ArgminAuglag = ctrlpp::argmin_solver<double, ctrlpp::argmin_auglag<>>;

auto bounded_settings() -> ctrlpp::argmin_settings<double>
{
    ctrlpp::argmin_settings<double> cfg{};
    cfg.max_time = 2.0;
    return cfg;
}

auto auglag_settings() -> ctrlpp::argmin_settings<double>
{
    ctrlpp::argmin_settings<double> cfg = bounded_settings();
    cfg.max_eval = 1000;
    return cfg;
}

// The variant list is spelled once; each pass over it supplies its own action.
// The two variants that regularly exhaust the time bound keep their single-shot
// quality record and stay out of the timed set, where nanobench's repetition
// would let them dominate the whole suite.
template <typename Action>
void for_each_variant(Action&& action)
{
    ctrlpp::nlopt_settings<double> nlopt_cfg{};
    nlopt_cfg.algorithm = ctrlpp::nlopt_algorithm::slsqp;
    const auto cfg = bounded_settings();
    action(NloptSolver{nlopt_cfg}, "nlopt", "slsqp", true, true);
    action(ArgminSlsqp{cfg}, "argmin", "slsqp", true, true);
    action(ArgminNwSqp{cfg}, "argmin", "nw_sqp", true, true);
    action(ArgminFilterSlsqp{cfg}, "argmin", "filter_slsqp", true, true);
    action(ArgminFilterNwSqp{cfg}, "argmin", "filter_nw_sqp", false, false);
    action(ArgminAuglag{auglag_settings()}, "argmin", "auglag", false, true);
}

// nanobench clears its accumulated results whenever the title changes, so a
// per-cell title would leave only the last cell in the rendered file; the cell
// rides in the row name instead. The family rides there too, because two arms
// both called "slsqp" would collide and a lookup by name would then hand one
// arm the other's figure.
auto row_label(char const* family, char const* algorithm, const std::string& system_name, std::size_t nx,
               int horizon) -> std::string
{
    return std::string{family} + "_" + algorithm + " " + system_name + " NX=" + std::to_string(nx) + " N="
         + std::to_string(horizon);
}

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_benchmark(const std::string& system_name, Dynamics dynamics, int horizon,
                   ankerl::nanobench::Bench& bench, std::ostream& quality_csv)
{
    const auto config = problems::make_nmpc_quadratic_config<NX, NU>(horizon);
    const auto x0 = problems::unit_first_axis_x0<NX>();
    const std::string warm_start = warm_start_label(bounded_settings().warm_start);

    std::vector<arms::arm_answer> answers;
    for_each_variant(
        [&](auto solver, char const* family, char const* algorithm, bool benched, bool)
        {
            const std::string label = row_label(family, algorithm, system_name, NX, horizon);
            auto probe = arms::probe_nmpc_arm<NX, NU>(dynamics, config, x0, std::move(solver), label);
            write_quality_csv_row(quality_csv, system_name, family, algorithm, warm_start,
                                  static_cast<int>(NX), horizon, probe.quality);
            if(benched)
                answers.push_back(probe.answer);
        });

    const double spread = arms::solution_spread(answers);
    bench.warmup((NX >= 8 && horizon >= 20) ? 5 : 50)
        .minEpochIterations((NX >= 8) ? ((horizon >= 20) ? 3 : 10) : 50);
    for_each_variant(
        [&](auto solver, char const* family, char const* algorithm, bool benched, bool)
        {
            const arms::arm_answer* answer =
                arms::answer_for_label(answers, row_label(family, algorithm, system_name, NX, horizon));
            if(!benched || answer == nullptr)
                return;
            auto controller = ctrlpp::bench::built_or_exit(
                ctrlpp::nmpc_dynamic<double, NX, NU, decltype(solver), Dynamics>::create(dynamics, config,
                                                                                         std::move(solver)),
                "arm");
            arms::emit_variant_rows(bench, *answer, spread,
                                    [&] { ankerl::nanobench::doNotOptimizeAway(controller.solve(x0)); });
        });
}

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_convergence(const std::string& system_name, Dynamics dynamics, int horizon, std::ostream& quality_csv)
{
    arms::write_convergence_rates<NX, NU>(system_name, dynamics,
                                          problems::make_nmpc_quadratic_config<NX, NU>(horizon), horizon,
                                          quality_csv, [](auto&& action) { for_each_variant(action); });
}

void run_size_sweep(ankerl::nanobench::Bench& bench, std::ostream& quality_csv)
{
    for(int horizon : {10, 20, 30})
    {
        run_benchmark<2, 1>("double_integrator", problems::double_integrator_2, horizon, bench, quality_csv);
        run_benchmark<4, 2>("double_integrator", problems::double_integrator_4, horizon, bench, quality_csv);
        run_benchmark<8, 4>("double_integrator", problems::double_integrator_8, horizon, bench, quality_csv);
    }
    for(int horizon : {10, 20})
        run_benchmark<2, 1>("pendulum", problems::pendulum_2, horizon, bench, quality_csv);
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("NMPC: sequential quadratic programming variants")
        .warmup(50)
        .minEpochIterations(50)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    std::ofstream timing_csv("bench_sqp_variants_timing.csv");
    std::ofstream quality_csv("bench_sqp_variants_quality.csv");
    write_quality_csv_header(quality_csv);

    run_size_sweep(bench, quality_csv);
    run_convergence<4, 2>("double_integrator", problems::double_integrator_4, 10, quality_csv);
    run_convergence<2, 1>("pendulum", problems::pendulum_2, 10, quality_csv);

    bench.render(ctrlpp::bench::csv_tpl, timing_csv);
}
