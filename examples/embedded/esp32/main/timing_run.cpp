#include "eigen_alloc_sentinel.h"

#include "telemetry.h"
#include "timing_run.h"
#include "cycle_record.h"
#include "cycle_counter.h"

#include "dsp_demo.h"
#include "estimation_demo.h"
#include "trajectory_demo.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_task_wdt.h"

#include <cstdio>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

namespace {

// A warm figure repeats its region as many times as the golden window has
// steps, so every warm worst-of-N shares N.
constexpr std::uint32_t kWarmRepeats = static_cast<std::uint32_t>(kGoldenSteps);

static_assert(kWarmRepeats <= kTimingCapacity);

// Static rather than on the control task's stack, which the figures share with
// the rest of the report.
cycle_record samples;

void report(const timing_region &region)
{
    const cycle_distribution d = samples.summarize();
    emit_spread_row(region, d);
    emit_interruption(region, d);
}

// The code and the constants it reads come from external flash through the
// cache, so the cold figure pays the misses a warm one may not; within the warm
// series the spread is interrupts landing in a sample and the other core
// competing for the same flash interface. The task yields between regions and
// never inside one.
template<class Region>
void time_region(const char *name, const char *cold_basis, std::uint32_t block, Region region)
{
    samples.clear();
    samples.add(region());
    report({name, "cold", cold_basis, block});

    samples.clear();
    for(std::uint32_t r = 0; r < kWarmRepeats; ++r)
        samples.add(region());
    report({name, "warm", "golden_window", block});
    vTaskDelay(1);
}

cycle_span time_design()
{
    const cycle_mark start = open_region();
    const auto demo        = control_loop_demo<float>::make();
    keep(demo);
    return close_region(start);
}

// Every block starts from the same copy, so the blocks repeat identical
// arithmetic and their spread is the hardware's and the scheduler's alone.
template<class Demo>
cycle_span time_block(const Demo &pristine)
{
    Demo demo = pristine;
    keep(demo);
    const cycle_mark start = open_region();
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
        keep_value(demo.step());
    keep(demo);
    return close_region(start);
}

template<template<class> class Demo>
void time_steps(const char *name)
{
    const auto demo = Demo<float>::make();
    if(!demo.has_value())
    {
        char body[kLineBytes];
        std::snprintf(body, sizeof(body), "timing region=%s verdict=REFUSED -- the family could not be built", name);
        send_line(stream_tag::verdict, body);
        return;
    }
    time_region(name, "first_after_reset", static_cast<std::uint32_t>(kGoldenSteps), [&demo] { return time_block(*demo); });
}

}

void stop_task_watchdog() noexcept
{
    esp_task_wdt_deinit();
}

counter_state run_timing()
{
    const counter_state state = prove_cycle_counter();
    emit_counter_canary(state);
    if(state != counter_state::live)
        return state;
    time_region("overhead", "first_after_proof", 1U, measure_overhead);
    time_region("control.design", "first_after_reset", 1U, time_design);
    time_steps<control_loop_demo>("control.step");
    time_steps<estimation_demo>("estimation.step");
    time_steps<dsp_demo>("dsp.step");
    time_steps<trajectory_demo>("trajectory.step");
    return state;
}

}
