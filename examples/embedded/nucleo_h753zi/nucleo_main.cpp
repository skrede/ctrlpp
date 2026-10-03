// Armed before Eigen is parsed, or the Eigen allocation sentinel is never compiled in.
#include "eigen_alloc_sentinel.h"

// ctrlpp NUCLEO-H753ZI on-device control-loop leg (bare superloop, double FPU).
//
// CMSIS startup runs SystemInit + C++ static ctors, then calls main(). This
// proves the allocation sensor is not blind, designs the infinite-horizon LQR
// gain once, runs the closed loop with the sensor armed around each step,
// streams the trajectory as CSV over USART3 -> ST-Link VCP, diffs the on-target
// double result against the independent host-double golden, and halts for the
// operator to read the console.

#include "build_id.h"
#include "alloc_sensor.h"
#include "sbrk_ceiling.h"

#include "golden_verdict.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <Eigen/Dense>

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp {
void usart3_console_init() noexcept;
}

namespace {

using demo_type = ctrlpp::control_loop_demo<double>;

ctrlpp::golden_bound golden_workspace = ctrlpp::make_golden_bound();

void halt() noexcept
{
    for(;;)
    {
    }
}

void report_identity()
{
    char build_id[ctrlpp::kBuildIdTextLength];
    ctrlpp::render_build_id(build_id);
    std::printf("[meta] build_id=%s\n", build_id);
    std::printf("[ctrlpp] NUCLEO-H753ZI bare-superloop control loop (double)\n");
}

void report_canary()
{
    const ctrlpp::canary_verdict verdict = ctrlpp::run_alloc_canary();
    std::printf("[canary] %s observed=%" PRIu32 " expected=%" PRIu32 " eigen_trips=%" PRIu32 "\n", verdict.live ? "PASS" : "FAIL", verdict.observed, verdict.expected,
                verdict.eigen_trips);
}

void report_heap()
{
    std::printf("[heap] setup_high_water=%" PRIu32 " reserve=%" PRIu32 " bytes\n", static_cast<std::uint32_t>(ctrlpp::heap_high_water_bytes()),
                static_cast<std::uint32_t>(ctrlpp::heap_reserve_bytes()));
}

void report_gain(const demo_type &demo)
{
    std::printf("gain    K = [%.9f, %.9f]\n", demo.K(0, 0), demo.K(0, 1));
    std::printf("host    K = [%.9f, %.9f]\n", ctrlpp::kHostK0, ctrlpp::kHostK1);
}

// The sensor is armed around step() alone, so the figure is the controller's
// and the telemetry formatting between steps is not charged to it.
void run_loop(demo_type &demo)
{
    std::printf("k,t_s,x0,x1,u\n");
    ctrlpp::alloc_sensor_reset();
    for(std::int32_t k = 0; k <= ctrlpp::kSteps; ++k)
    {
        ctrlpp::alloc_sensor_arm();
        const double u = demo.step();
        ctrlpp::alloc_sensor_disarm();
        std::printf("%" PRId32 ",%.4f,%.6f,%.6f,%.6f\n", k, k * ctrlpp::kDt, demo.x[0], demo.x[1], u);
    }
}

void report_allocations()
{
    const std::uint32_t observed = ctrlpp::alloc_sensor_observed();
    const std::uint32_t trips    = ctrlpp::alloc_sensor_eigen_trips();
    const std::int32_t steps     = ctrlpp::kSteps + 1;
    std::printf("[alloc] steady_state allocations=%" PRIu32 " steps=%" PRId32 " per_step=%.2f eigen_sentinel=%s trips=%" PRIu32 "\n", observed, steps,
                static_cast<double>(observed) / steps, trips == 0 ? "clean" : "TRIPPED", trips);
}

void report_golden(const demo_type &demo)
{
    const ctrlpp::golden_verdict verdict = ctrlpp::judge_golden(demo, golden_workspace);
    std::printf("gain check %s: dev = [%.3e, %.3e], tol = %.3e\n", verdict.gain_pass ? "PASS" : "FAIL", verdict.gain_departure[0], verdict.gain_departure[1],
                verdict.gain_tolerance);
    std::printf("settling: final |x| = %.6e (host = %.6e, err = %.3e, tol = %.3e)\n", demo.x.norm(), ctrlpp::kHostFinalNorm, verdict.norm_departure,
                verdict.norm_tolerance);
    std::printf("golden diff %s (gain %s, final norm %s)\n", verdict.pass ? "PASS" : "FAIL", verdict.gain_pass ? "PASS" : "FAIL", verdict.norm_pass ? "PASS" : "FAIL");
}

}

int main()
{
    ctrlpp::usart3_console_init();
    report_identity();
    report_canary();

    auto demo = demo_type::make();
    if(!demo.has_value())
    {
        std::printf("[ctrlpp] lqr_gain REFUSED the plant on device -- %s\n", ctrlpp::describe(demo.error()));
        halt();
    }

    report_heap();
    report_gain(*demo);
    run_loop(*demo);
    report_allocations();
    report_golden(*demo);
    halt();
}
