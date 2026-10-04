#include "eigen_alloc_sentinel.h"

#include "alloc_sensor.h"
#include "sbrk_ceiling.h"
#include "family_report.h"

#include "golden_reference.h"

#if CTRLPP_MCU_PREDICTIVE
    #include "predictive_demo.h"
    #include "predictive_verdict.h"
#endif

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp {

window_figures read_window() noexcept
{
    return {alloc_sensor_observed(), alloc_sensor_eigen_trips()};
}

void report_allocations(const window_figures &figures)
{
    const std::uint32_t steps = static_cast<std::uint32_t>(kGoldenSteps);
    std::printf("[alloc] steady_state allocations=%" PRIu32 " steps=%" PRIu32 " per_step=%.2f eigen_sentinel=%s trips=%" PRIu32 "\n", figures.allocations, steps,
                static_cast<double>(figures.allocations) / steps, figures.eigen_trips == 0 ? "clean" : "TRIPPED", figures.eigen_trips);
}

void report_family(const family_result &result, const window_figures &figures)
{
    const std::uint32_t steps = static_cast<std::uint32_t>(kGoldenSteps);
    std::printf("[family] %s value=%.17g golden=%.17g departure=%.3e bound=%.3e verdict=%s allocations=%" PRIu32 " steps=%" PRIu32 " per_step=%.2f eigen_trips=%" PRIu32 "\n",
                result.name, result.value, result.golden, result.departure, result.bound, result.pass ? "PASS" : "FAIL", figures.allocations, steps,
                static_cast<double>(figures.allocations) / steps, figures.eigen_trips);
}

#if CTRLPP_MCU_PREDICTIVE

namespace {

// The controller's constructor materializes a second solver of about 41 KB and
// moves it in; building the instance in a call that returns first pops that
// copy off the stack before the solves need the room.
[[gnu::noinline]] predictive_demo<double> &predictive_instance()
{
    static predictive_demo<double> demo;
    return demo;
}

predictive_record predictive_run{};
predictive_workspace predictive_bound = make_predictive_workspace();

void report_predictive(const predictive_verdict &verdict)
{
    std::printf("[predictive] stationarity_stops=%" PRIu32 " of %" PRIu32 " worst_gate=%.3e run_steps=%" PRIu32 " armed_steps=%" PRIu32 " premise=%s gate=%s\n",
                predictive_run.stationarity_stops, static_cast<std::uint32_t>(kPredictiveRunSteps), verdict.worst_gate, static_cast<std::uint32_t>(kPredictiveRunSteps),
                static_cast<std::uint32_t>(kGoldenSteps), verdict.premise ? "PASS" : "FAIL", verdict.gate ? "PASS" : "FAIL");
}

}

// The armed steps allocate nothing, so the high-water after the run is the
// setup's, construction and the first solves included.
void drive_predictive_family()
{
    predictive_demo<double> &demo = predictive_instance();
    alloc_sensor_reset();
    const auto ran               = run_predictive(demo, predictive_run, alloc_sensor_arm, alloc_sensor_disarm);
    const window_figures figures = read_window();
    std::printf("[heap] predictive_setup_high_water=%" PRIu32 " reserve=%" PRIu32 " bytes\n", static_cast<std::uint32_t>(heap_high_water_bytes()),
                static_cast<std::uint32_t>(heap_reserve_bytes()));
    if(!ran.has_value())
    {
        std::printf("[family] predictive verdict=REFUSED -- %s\n", describe(ran.error()));
        return;
    }
    const predictive_verdict verdict = judge_predictive(predictive_run, predictive_bound);
    report_predictive(verdict);
    report_family({"predictive", verdict.cost, kHostPredictiveCost, verdict.departure, verdict.bound, verdict.pass}, figures);
}

#else

void drive_predictive_family()
{
    std::printf("[family] predictive verdict=ABSENT -- the image was built without the predictive backend\n");
}

#endif

}
