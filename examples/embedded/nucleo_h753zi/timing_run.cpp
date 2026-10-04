#include "eigen_alloc_sentinel.h"

#include "timing_run.h"
#include "cycle_record.h"
#include "cycle_counter.h"

#include "dsp_demo.h"
#include "estimation_demo.h"
#include "trajectory_demo.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cinttypes>
#include <algorithm>

namespace ctrlpp {

namespace {

// A warm figure repeats its region as many times as the golden window has
// steps, so every warm worst-of-N, the predictive solves' included, shares N.
constexpr std::uint32_t kWarmRepeats = static_cast<std::uint32_t>(kGoldenSteps);

static_assert(kWarmRepeats <= kTimingCapacity);

template<class Region>
void time_region(const char *name, const char *cold_basis, std::uint32_t block, Region region)
{
    cycle_record cold;
    cold.add(region());
    report_timing({name, "cold", cold_basis, block}, cold.summarize());

    cycle_record warm;
    for(std::uint32_t r = 0; r < kWarmRepeats; ++r)
        warm.add(region());
    report_timing({name, "warm", "golden_window", block}, warm.summarize());
}

cycle_span time_design()
{
    const std::uint32_t start = cycle_sample();
    const auto demo           = control_loop_demo<double>::make();
    keep(demo);
    return cycles_between(start, cycle_sample());
}

// Every block starts from the same copy, so the blocks repeat identical
// arithmetic and their spread is the hardware's alone.
template<class Demo>
cycle_span time_block(const Demo &pristine)
{
    Demo demo = pristine;
    keep(demo);
    const std::uint32_t start = cycle_sample();
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
        keep_value(demo.step());
    keep(demo);
    return cycles_between(start, cycle_sample());
}

template<template<class> class Demo>
void time_steps(const char *name)
{
    const auto demo = Demo<double>::make();
    if(!demo.has_value())
    {
        std::printf("[timing] region=%s verdict=REFUSED -- the family could not be built\n", name);
        return;
    }
    time_region(name, "first_after_reset", static_cast<std::uint32_t>(kGoldenSteps), [&demo] { return time_block(*demo); });
}

}

counter_state run_timing()
{
    const counter_state state = enable_cycle_counter();
    if(state != counter_state::live)
    {
        std::printf("[timing] UNAVAILABLE -- %s; no timing figure is printed\n", describe(state));
        return state;
    }
    time_region("overhead", "first_after_proof", 1U, measure_overhead);
    time_region("control.design", "first_after_reset", 1U, time_design);
    time_steps<control_loop_demo>("control.step");
    time_steps<estimation_demo>("estimation.step");
    time_steps<dsp_demo>("dsp.step");
    time_steps<trajectory_demo>("trajectory.step");
    return state;
}

solve_timer::solve_timer()
        : spans_{}
        , iterations_{}
        , start_(0)
{
}

void solve_timer::enter(std::size_t) noexcept
{
    start_ = cycle_sample();
}

void solve_timer::leave(std::size_t k, std::int32_t iterations) noexcept
{
    spans_[k]      = cycles_between(start_, cycle_sample());
    iterations_[k] = iterations;
}

void solve_timer::report() const
{
    if(prove_cycle_counter() != counter_state::live)
        return;
    cycle_record cold;
    cold.add(spans_[0]);
    report_timing({"predictive.step", "cold", "first_after_reset", 1U}, cold.summarize());

    cycle_record warm;
    for(std::size_t k = static_cast<std::size_t>(kPredictiveWarmupSolves); k < kPredictiveRunSteps; ++k)
        warm.add(spans_[k]);
    const cycle_distribution window = warm.summarize();
    report_timing({"predictive.step", "warm", "golden_window", 1U}, window);
    report_work(window.max_at);
}

// The solves differ in work as well as in fetch state, so the iteration counts
// are printed beside the figures they explain.
void solve_timer::report_work(std::uint32_t slowest_at) const
{
    const auto first  = iterations_.begin() + kPredictiveWarmupSolves;
    const auto bounds = std::minmax_element(first, iterations_.end());
    std::printf("[work] region=predictive.step regime=cold iterations=%" PRId32 "\n", iterations_[0]);
    std::printf("[work] region=predictive.step regime=warm n=%" PRIu32 " iterations_min=%" PRId32 " iterations_max=%" PRId32 " max_at_iterations=%" PRId32 "\n",
                static_cast<std::uint32_t>(iterations_.end() - first), *bounds.first, *bounds.second, *(first + slowest_at));
}

}
