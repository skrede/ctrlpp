// Competitive benchmark: ctrlpp::care vs ct::optcon::CARE
// Problem: continuous-time Riccati solve for damped chain-of-integrators systems,
// size-swept NX in {2, 4, 6, 8, 12, 16, 20, 24, 30}. NU scales with NX.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ct_care_arm.h"

#include <fstream>

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("CARE: ctrlpp vs ct_optcon (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    using ctrlpp::bench::run_care_sweep;
    run_care_sweep<2, 1>(bench,  "ctrlpp::care NX=2",  "ct::optcon::CARE NX=2");
    run_care_sweep<4, 2>(bench,  "ctrlpp::care NX=4",  "ct::optcon::CARE NX=4");
    run_care_sweep<6, 2>(bench,  "ctrlpp::care NX=6",  "ct::optcon::CARE NX=6");
    run_care_sweep<8, 2>(bench,  "ctrlpp::care NX=8",  "ct::optcon::CARE NX=8");
    run_care_sweep<12, 3>(bench, "ctrlpp::care NX=12", "ct::optcon::CARE NX=12");
    run_care_sweep<16, 4>(bench, "ctrlpp::care NX=16", "ct::optcon::CARE NX=16");
    run_care_sweep<20, 5>(bench, "ctrlpp::care NX=20", "ct::optcon::CARE NX=20");
    run_care_sweep<24, 6>(bench, "ctrlpp::care NX=24", "ct::optcon::CARE NX=24");
    run_care_sweep<30, 6>(bench, "ctrlpp::care NX=30", "ct::optcon::CARE NX=30");

    std::ofstream csv("bench_care_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
