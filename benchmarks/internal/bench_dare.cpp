// Discrete Riccati benchmark, swept over state dimension.
//
// The acceptance check is the independent variable. Four quantities are timed
// per size; the first three run on identical operands:
//
//   solve+accept -- the whole `ctrlpp::dare` entry point
//   accept       -- `verify_dare_solution` alone, on the solve's own answer
//   fwd-error    -- `estimate_riccati_forward_error` alone, on pre-formed
//                   closed-loop and residual operands
//   solve+accept scaled
//                -- the same entry point on weights whose scale is s rather
//                   than one, timed at every rung of a geometric ladder in s
//
// The solve without the check is the difference of the first two, which is the
// denominator the published operation counts use: those are stated against a
// solve with no acceptance arithmetic in it. Timing the round's start against
// its end would confound the construction change, the removed factorization and
// the operand preparation into one number that attributes to none of them.
//
// The corpus is the discrete damped chain: a forward-Euler step of the same
// continuous chain the CARE bakeoff sweeps, so the two harnesses describe the
// same family. It is built at two weight scales. Q = I and R = 0.1 I put the
// weight scale at exactly one, so equilibration is inactive and the solve runs
// the acceptance check ONCE. Multiplying both weights by s puts the weight
// scale at exactly s, and the entry point then runs the whole equilibrated
// path: the acceptance check at the equilibrated scale, the rescale of the
// answer, a definiteness factorization of the rescaled answer, the acceptance
// check again at the caller's scale, and a gain agreement between the two
// scales. That is the fourth quantity above.
//
// THE EQUILIBRATED PATH'S COST IS THE DIFFERENCE between the scaled row and the
// scale-one row at the same size. It is not timed separately and it is not
// asserted anywhere, so a reader recomputes it from two published rows rather
// than taking one number on trust.
//
// `--check` runs the corpus confirmation instead of the timing set: it prints,
// for every size and rung, the weight scale the helper computes and whether the
// entry point's equilibrated gate is therefore open. A scaled corpus that
// failed to open that gate would produce a second row equal to the first and an
// overhead of approximately zero, and nothing in a timing would show it.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "dare_rows.h"
#include "dare_corpus.h"
#include "dare_corpus_check.h"

#include <cstdio>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string_view>

namespace
{

using ctrlpp::bench::for_each_size;
using ctrlpp::bench::run_size_sweep;

/// AN ARGUMENT THIS PROGRAM DOES NOT RECOGNIZE MUST NOT FALL THROUGH TO THE
/// TIMING SET. The timed rows are comparable to the published figures only when
/// the machine was quiet, and an argument typed by someone expecting a usage
/// message is exactly the invocation that would otherwise produce a timing
/// nobody intended, on a machine nobody had reserved. The measurement runs on no
/// arguments and on nothing else.
int usage(char const* unrecognized)
{
    if (unrecognized != nullptr)
        std::fprintf(stderr, "bench_dare: unrecognized argument '%s'\n", unrecognized);
    std::fprintf(unrecognized != nullptr ? stderr : stdout,
                 "usage: bench_dare           run the timing set (reserved machine only)\n"
                 "       bench_dare --check   confirm the scaled corpus, timing nothing\n"
                 "       bench_dare --smoke   run the timing set at one iteration, measuring nothing\n"
                 "       bench_dare --help    this message\n");
    return unrecognized != nullptr ? 2 : 0;
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
    bool check_only = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string_view arg(argv[i]);
        if (arg == "--check")
            check_only = true;
        else if (arg == "--smoke")
            continue;
        else if (arg == "--help" || arg == "-h")
            return usage(nullptr);
        else
            return usage(argv[i]);
    }

    if (check_only)
        return ctrlpp::bench::run_corpus_check();

    check_perf_event_paranoid();

    // 51 epochs give an odd-sized sample, so the reported median is an observed
    // measurement rather than an interpolation, and the CSV's error column is
    // the median absolute percent error across those same 51.
    ankerl::nanobench::Bench bench;
    bench.title("DARE acceptance overhead (size sweep)")
        .warmup(10)
        .epochs(51)
        .minEpochTime(std::chrono::milliseconds(1))
        .performanceCounters(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    for_each_size([&]<std::size_t NX, std::size_t NU>() { run_size_sweep<NX, NU>(bench); });

    std::ofstream csv("bench_dare.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
