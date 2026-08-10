// Competitive benchmark: ctrlpp::care vs ct::optcon::CARE at NX=30, NU=6.
// Single-size profile bench used as a perf target under -fno-omit-frame-pointer.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "comparison/ct/ct_care_arm.h"

#include <fstream>

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("CARE profile bench NX=30")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    ctrlpp::bench::run_care_sweep<30, 6>(bench, "ctrlpp::care NX=30", "ct::optcon::CARE NX=30");

    std::ofstream csv("bench_care_vs_ct_nx30.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
