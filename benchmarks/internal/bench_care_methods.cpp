// Bakeoff benchmark: race three ctrlpp::care method tags
// (schur, sign_function, balanced_schur) against ct::optcon::CARE on
// identical damped-chain inputs, size-swept NX in {2, 4, 6, 8, 12, 16,
// 20, 24, 30}.
//
// Primary metric: median instructions (nanobench performanceCounters).
// Secondary: wall-clock. Winner selection is mechanical; see tools/bakeoff_winner.py.
//
// The sweep also carries ONE acceptance-check row per size. Every method tag
// reaches the same postcondition, against the same Hamiltonian, so its cost is
// a single measurement -- but each tag's solve is a different denominator, so
// the one row lands in three overhead statements. The published continuous
// count is stated against the solve alone, so the timed denominator is the tag's
// row minus the acceptance row rather than the tag's row itself.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "comparison/ct/ct_care_arm.h"

#include "ctrlpp/detail/care_postconditions.h"

#include <cstdio>
#include <chrono>
#include <string>
#include <fstream>
#include <cstddef>
#include <iostream>

namespace
{

using ctrlpp::bench::care_residual_metric;
using ctrlpp::bench::riccati_plant;

// Four arms admit no pairwise agreement figure, so what the arm rows carry is
// the spread of the set: how far the widest-separated pair of the four
// solutions stand apart, one number shared by all four rows exactly as a
// two-arm deviation is shared by both of its rows.
constexpr char const* spread_metric = "max abs entrywise spread of the four arms' Riccati solutions P";

template <std::size_t NX>
using care_solution = Eigen::Matrix<double, int(NX), int(NX)>;

template <std::size_t NX, std::size_t NU>
struct method_arms
{
    care_solution<NX> schur;
    care_solution<NX> sign;
    care_solution<NX> balanced;
    care_solution<NX> competitor;
};

template <std::size_t NX, std::size_t NU>
method_arms<NX, NU> solve_all(const riccati_plant<NX, NU>& plant, ctrlpp::bench::ct_care_arm<NX, NU>& ct_arm)
{
    using ctrlpp::bench::built_or_exit;
    return {built_or_exit(ctrlpp::care<double, NX, NU, ctrlpp::detail::schur_care_method>(
                              plant.A, plant.B, plant.Q, plant.R), "ctrlpp::care[schur]").P,
            built_or_exit(ctrlpp::care<double, NX, NU, ctrlpp::detail::sign_function_care_method>(
                              plant.A, plant.B, plant.Q, plant.R), "ctrlpp::care[sign]").P,
            built_or_exit(ctrlpp::care<double, NX, NU, ctrlpp::detail::balanced_schur_care_method>(
                              plant.A, plant.B, plant.Q, plant.R), "ctrlpp::care[balanced]").P,
            ct_arm.solve()};
}

template <std::size_t NX, std::size_t NU>
double entrywise_spread(const method_arms<NX, NU>& arms)
{
    const care_solution<NX> high =
        arms.schur.cwiseMax(arms.sign).cwiseMax(arms.balanced).cwiseMax(arms.competitor);
    const care_solution<NX> low =
        arms.schur.cwiseMin(arms.sign).cwiseMin(arms.balanced).cwiseMin(arms.competitor);
    return (high - low).cwiseAbs().maxCoeff();
}

template <std::size_t NX, std::size_t NU, typename Op>
void emit_arm(ankerl::nanobench::Bench& bench, const std::string& label, const riccati_plant<NX, NU>& plant,
              const care_solution<NX>& P, double spread, Op&& op)
{
    ctrlpp::bench::run_with_accuracy(bench, spread_metric, label, spread, op);
    ctrlpp::bench::run_own_criterion_row(bench, care_residual_metric, label.c_str(),
                                         ctrlpp::bench::riccati_relative_residual<NX, NU>(plant, P),
                                         std::forward<Op>(op));
}

template <std::size_t NX, std::size_t NU, typename Method>
auto tag_solver(const riccati_plant<NX, NU>& plant)
{
    return [&plant] {
        auto r = ctrlpp::care<double, NX, NU, Method>(plant.A, plant.B, plant.Q, plant.R);
        ankerl::nanobench::doNotOptimizeAway(r);
    };
}

// The acceptance check, timed on the operands the solver hands it: the
// Hamiltonian the caller's problem defines and the matrix the extraction
// produced. Q = I and R = 0.1 I put the weight scale at exactly one, so the
// equilibrated Hamiltonian and the caller's are the same object here and this
// row measures what every tag actually pays. It answers a verdict rather than a
// Riccati equation, so it has no residual of its own.
template <std::size_t NX, std::size_t NU>
void emit_accept_row(ankerl::nanobench::Bench& bench, const riccati_plant<NX, NU>& plant,
                     const care_solution<NX>& P)
{
    auto H = ctrlpp::detail::build_care_hamiltonian<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R);
    if(!H)
    {
        std::cerr << "SKIP acceptance row NX=" << NX << ": operands unavailable\n";
        return;
    }
    const Eigen::Matrix<double, 2 * int(NX), 2 * int(NX)> H_accept = *H;
    if(!ctrlpp::detail::care_solution_satisfies_postconditions<double, NX>(H_accept, P))
    {
        std::cerr << "SKIP acceptance row NX=" << NX << ": postcondition declined\n";
        return;
    }

    char buf[64];
    std::snprintf(buf, sizeof(buf), "ctrlpp::care[accept] NX=%zu", NX);
    ctrlpp::bench::run_single_implementation_row(bench, buf, [&] {
        auto ok = ctrlpp::detail::care_solution_satisfies_postconditions<double, NX>(H_accept, P);
        ankerl::nanobench::doNotOptimizeAway(ok);
    });
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench)
{
    const riccati_plant<NX, NU> plant = ctrlpp::bench::build_damped_chain<NX, NU>();
    ctrlpp::bench::ct_care_arm<NX, NU> ct_arm{plant};
    const method_arms<NX, NU> arms = solve_all<NX, NU>(plant, ct_arm);
    const double spread = entrywise_spread<NX, NU>(arms);
    const std::string size = " NX=" + std::to_string(NX);

    emit_arm<NX, NU>(bench, "ctrlpp::care[schur]" + size, plant, arms.schur, spread,
                     tag_solver<NX, NU, ctrlpp::detail::schur_care_method>(plant));
    emit_arm<NX, NU>(bench, "ctrlpp::care[sign]" + size, plant, arms.sign, spread,
                     tag_solver<NX, NU, ctrlpp::detail::sign_function_care_method>(plant));
    emit_arm<NX, NU>(bench, "ctrlpp::care[balanced]" + size, plant, arms.balanced, spread,
                     tag_solver<NX, NU, ctrlpp::detail::balanced_schur_care_method>(plant));
    emit_arm<NX, NU>(bench, "ct::optcon::CARE" + size, plant, arms.competitor, spread, [&ct_arm] {
        auto P = ct_arm.solve();
        ankerl::nanobench::doNotOptimizeAway(P);
    });
    emit_accept_row<NX, NU>(bench, plant, arms.schur);
}

void check_perf_event_paranoid()
{
    std::ifstream paranoid_file("/proc/sys/kernel/perf_event_paranoid");
    int paranoid_value = 3;
    paranoid_file >> paranoid_value;
    if (paranoid_value > 2)
        std::cerr << "WARN: perf_event_paranoid=" << paranoid_value
                  << "; nanobench hardware counters may be zero. "
                  << "Run:  sudo sysctl kernel.perf_event_paranoid=2\n";
}

}

int main(int argc, char** argv)
{
    check_perf_event_paranoid();

    // 51 epochs give an odd-sized sample, so the reported median is an observed
    // measurement rather than an interpolation, and the CSV's error column is
    // the median absolute percent error across those same 51. Iterations per
    // epoch are set by a 1 ms floor rather than a fixed count, so a row at
    // NX = 30 does not cost three orders of magnitude more wall time than a row
    // at NX = 2 to reach the same timer resolution.
    ankerl::nanobench::Bench bench;
    bench.title("CARE methods bakeoff (size sweep)")
        .warmup(50)
        .epochs(51)
        .minEpochTime(std::chrono::milliseconds(1))
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    run_size_sweep<2,  1>(bench);
    run_size_sweep<4,  2>(bench);
    run_size_sweep<6,  2>(bench);
    run_size_sweep<8,  2>(bench);
    run_size_sweep<12, 3>(bench);
    run_size_sweep<16, 4>(bench);
    run_size_sweep<20, 5>(bench);
    run_size_sweep<24, 6>(bench);
    run_size_sweep<30, 6>(bench);

    std::ofstream csv("bench_care_methods.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
