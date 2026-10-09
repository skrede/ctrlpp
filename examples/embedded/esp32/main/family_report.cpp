#include "eigen_alloc_sentinel.h"

#include "telemetry.h"
#include "alloc_sensor.h"
#include "family_report.h"

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp {

namespace {

// The free-heap view rides beside the counter's figure rather than replacing
// it: two independent numbers that disagree are a finding, and neither is
// picked over the other.
void append_figures(char *body, std::size_t size, const window_figures &figures)
{
    const std::uint32_t steps = static_cast<std::uint32_t>(kGoldenSteps);
    std::snprintf(body, size, " allocations=%" PRIu32 " steps=%" PRIu32 " per_step=%.2f eigen_trips=%" PRIu32 " heap_retained=%" PRId32 " heap_dip=%" PRIu32, figures.allocations, steps,
                  static_cast<double>(figures.allocations) / steps, figures.eigen_trips, figures.heap.retained, figures.heap.dip);
}

}

window_figures read_window() noexcept
{
    if constexpr(kAllocSensorPresent)
        return {true, alloc_sensor_observed(), alloc_sensor_eigen_trips(), alloc_sensor_heap()};
    else
        return {false, 0, 0, {0, 0}};
}

void report_allocations(const window_figures &figures)
{
    if(!figures.measured)
    {
        send_line(stream_tag::verdict, "alloc steady_state allocations=unmeasured -- this image carries no allocation sensor");
        return;
    }
    const std::uint32_t steps = static_cast<std::uint32_t>(kGoldenSteps);
    char body[kLineBytes];
    std::snprintf(body, sizeof(body),
                  "alloc steady_state allocations=%" PRIu32 " steps=%" PRIu32 " per_step=%.2f eigen_sentinel=%s trips=%" PRIu32 " heap_retained=%" PRId32 " heap_dip=%" PRIu32,
                  figures.allocations, steps, static_cast<double>(figures.allocations) / steps, figures.eigen_trips == 0 ? "clean" : "TRIPPED", figures.eigen_trips,
                  figures.heap.retained, figures.heap.dip);
    send_line(stream_tag::verdict, body);
}

void report_family(const family_result &result, const window_figures &figures)
{
    char body[kLineBytes];
    const int used       = std::snprintf(body, sizeof(body), "family %s value=%.17g golden=%.17g departure=%.3e bound=%.3e verdict=%s", result.name, result.value, result.golden,
                                         result.departure, result.bound, result.pass ? "PASS" : "FAIL");
    const std::size_t at = used > 0 && static_cast<std::size_t>(used) < sizeof(body) ? static_cast<std::size_t>(used) : sizeof(body) - 1;
    if(figures.measured)
        append_figures(body + at, sizeof(body) - at, figures);
    else
        std::snprintf(body + at, sizeof(body) - at, " allocations=unmeasured");
    send_line(stream_tag::verdict, body);
}

void report_refusal(const char *name, const char *reason)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "family %s verdict=REFUSED -- %s", name, reason);
    send_line(stream_tag::verdict, body);
}

void report_absent_predictive()
{
    send_line(stream_tag::verdict, "family predictive verdict=ABSENT -- the image was built without the predictive backend");
}

void report_golden(const control_loop_demo<float> &demo, const golden_verdict &verdict, const window_figures &figures)
{
    const double final_norm = static_cast<double>(demo.x.norm());
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "gain check %s: dev = [%.3e, %.3e], tol = %.3e", verdict.gain_pass ? "PASS" : "FAIL", verdict.gain_departure[0], verdict.gain_departure[1],
                  verdict.gain_tolerance);
    send_line(stream_tag::verdict, body);
    std::snprintf(body, sizeof(body), "settling: final |x| = %.6e (host = %.6e, err = %.3e, tol = %.3e)", final_norm, kHostFinalNorm, verdict.norm_departure, verdict.norm_tolerance);
    send_line(stream_tag::verdict, body);
    report_family({"control", final_norm, kHostFinalNorm, verdict.norm_departure, verdict.norm_tolerance, verdict.pass}, figures);
    std::snprintf(body, sizeof(body), "golden diff %s (gain %s, final norm %s)", verdict.pass ? "PASS" : "FAIL", verdict.gain_pass ? "PASS" : "FAIL",
                  verdict.norm_pass ? "PASS" : "FAIL");
    send_line(stream_tag::verdict, body);
}

}
