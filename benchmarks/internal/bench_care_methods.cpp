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

#include "ctrlpp/detail/care_postconditions.h"
#include "ctrlpp/control/care.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>

#include <Eigen/Dense>

#include <cstdio>
#include <chrono>
#include <fstream>
#include <cstddef>
#include <iostream>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

template <std::size_t NX, std::size_t NU>
auto build_damped_chain()
{
    // Continuous-time damped chain: A has -0.5 on the diagonal and 1.0 on the superdiagonal.
    // Spectrum is Re(lambda) < 0 so CARE is well-defined.
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Zero();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) = -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = 1.0;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = 1.0;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();

    return std::tuple{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench)
{
    auto [A, B, Q, R] = build_damped_chain<NX, NU>();

    ct::optcon::CARE<NX, NU> ct_care;
    typename ct::optcon::CARE<NX, NU>::state_matrix_t        A_ct = A;
    typename ct::optcon::CARE<NX, NU>::control_gain_matrix_t B_ct = B;
    typename ct::optcon::CARE<NX, NU>::state_matrix_t        Q_ct = Q;
    typename ct::optcon::CARE<NX, NU>::control_matrix_t      R_ct = R;

    // Warm-ups (outside measurement window).
    auto w_schur = ctrlpp::care<double, NX, NU, ctrlpp::detail::schur_care_method>(A, B, Q, R);
    auto w_sign  = ctrlpp::care<double, NX, NU, ctrlpp::detail::sign_function_care_method>(A, B, Q, R);
    auto w_bal   = ctrlpp::care<double, NX, NU, ctrlpp::detail::balanced_schur_care_method>(A, B, Q, R);
    auto w_ct    = ct_care.computeSteadyStateRiccatiMatrix(Q_ct, R_ct, A_ct, B_ct);
    ankerl::nanobench::doNotOptimizeAway(w_schur);
    ankerl::nanobench::doNotOptimizeAway(w_sign);
    ankerl::nanobench::doNotOptimizeAway(w_bal);
    ankerl::nanobench::doNotOptimizeAway(w_ct);

    char buf[64];

    std::snprintf(buf, sizeof(buf), "ctrlpp::care[schur] NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto r = ctrlpp::care<double, NX, NU, ctrlpp::detail::schur_care_method>(A, B, Q, R);
        ankerl::nanobench::doNotOptimizeAway(r);
    });

    std::snprintf(buf, sizeof(buf), "ctrlpp::care[sign] NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto r = ctrlpp::care<double, NX, NU, ctrlpp::detail::sign_function_care_method>(A, B, Q, R);
        ankerl::nanobench::doNotOptimizeAway(r);
    });

    std::snprintf(buf, sizeof(buf), "ctrlpp::care[balanced] NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto r = ctrlpp::care<double, NX, NU, ctrlpp::detail::balanced_schur_care_method>(A, B, Q, R);
        ankerl::nanobench::doNotOptimizeAway(r);
    });

    std::snprintf(buf, sizeof(buf), "ct::optcon::CARE NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto P_ct = ct_care.computeSteadyStateRiccatiMatrix(Q_ct, R_ct, A_ct, B_ct);
        ankerl::nanobench::doNotOptimizeAway(P_ct);
    });

    // The acceptance check, timed on the operands the solver hands it: the
    // Hamiltonian the caller's problem defines and the matrix the extraction
    // produced. Q = I and R = 0.1 I put the weight scale at exactly one, so the
    // equilibrated Hamiltonian and the caller's are the same object here and
    // this row measures what every tag actually pays.
    auto H = ctrlpp::detail::build_care_hamiltonian<double, NX, NU>(A, B, Q, R);
    if (!H || !w_schur)
    {
        std::cerr << "SKIP acceptance row NX=" << NX << ": operands unavailable\n";
        return;
    }
    const Eigen::Matrix<double, 2 * int(NX), 2 * int(NX)> H_accept = *H;
    const Eigen::Matrix<double, int(NX), int(NX)>         P_accept = w_schur->P;
    if (!ctrlpp::detail::care_solution_satisfies_postconditions<double, NX>(H_accept, P_accept))
    {
        // A check that declines exits early and would time a fraction of the
        // work the accepted path does. Refusing to report it is the point.
        std::cerr << "SKIP acceptance row NX=" << NX << ": postcondition declined\n";
        return;
    }

    std::snprintf(buf, sizeof(buf), "ctrlpp::care[accept] NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto ok = ctrlpp::detail::care_solution_satisfies_postconditions<double, NX>(H_accept, P_accept);
        ankerl::nanobench::doNotOptimizeAway(ok);
    });
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

int main()
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
    bench.render(comma_csv_tpl, csv);
}
