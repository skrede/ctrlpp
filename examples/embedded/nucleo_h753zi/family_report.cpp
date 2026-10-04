#include "eigen_alloc_sentinel.h"

#include "alloc_sensor.h"
#include "family_report.h"

#include "golden_reference.h"

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

}
