// Armed before Eigen is parsed, or the Eigen allocation sentinel is never compiled in.
#include "eigen_alloc_sentinel.h"

// ctrlpp ESP32 on-device leg (FreeRTOS, single-precision FPU).
//
// Drives the shared float control kernel from a fixed-cadence FreeRTOS task on
// real Xtensa silicon: prints the reset that started the run and the posture it
// is taken under, then either proves the allocation sensor is not blind (the
// evidence image) or proves the cycle counter live and times every family cold
// and warm (the timing image). It then designs the infinite-horizon LQR gain
// once, runs the closed loop at 50 Hz with the sensor armed around each step
// while streaming the trajectory's exact rows over UART2 (telemetry adapter,
// ESP TX GPIO5 -> USB RX), runs the estimation, signal-processing and trajectory
// families each inside its own armed window, and diffs every family's float
// result against its host-double golden.

#include "posture.h"
#include "telemetry.h"
#include "timing_run.h"
#include "stack_probe.h"
#include "alloc_sensor.h"
#include "reset_reason.h"
#include "family_report.h"

#include "dsp_demo.h"
#include "golden_verdict.h"
#include "estimation_demo.h"
#include "trajectory_demo.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_log.h"

#include <cstdint>

namespace {

using demo_type = ctrlpp::control_loop_demo<float>;

constexpr const char *TAG = "ctrlpp_esp32";

ctrlpp::golden_bound golden_workspace = ctrlpp::make_golden_bound();

// Runs before anything else executes a timed region, so the timing image's
// cold figures are first executions after reset.
void prove_instrument()
{
    if constexpr(ctrlpp::kAllocSensorPresent)
        ctrlpp::emit_allocation_canary(ctrlpp::run_alloc_canary());
    else
        ctrlpp::run_timing();
}

void report_opening()
{
    ctrlpp::emit_identity();
    ctrlpp::emit_posture(ctrlpp::read_reset_cause());
    ctrlpp::emit_header(ctrlpp::claim_class::spread);
    prove_instrument();
}

// The sensor is armed around step() alone, so the figure is the controller's
// and the telemetry and the cadence wait between steps are not charged to it.
void run_loop(demo_type &demo)
{
    ctrlpp::emit_header(ctrlpp::claim_class::exact);
    ctrlpp::alloc_sensor_reset();
    TickType_t next         = xTaskGetTickCount();
    const TickType_t period = pdMS_TO_TICKS(20);
    for(std::int32_t k = 0; k <= ctrlpp::kSteps; ++k)
    {
        ctrlpp::alloc_sensor_arm();
        const float u = demo.step();
        ctrlpp::alloc_sensor_disarm();
        ctrlpp::emit_exact_row(k, {static_cast<float>(k * ctrlpp::kDt), demo.x[0], demo.x[1], u});
        vTaskDelayUntil(&next, period);
    }
}

void run_families(const ctrlpp::window_figures &control)
{
    ctrlpp::report_allocations(control);
    ctrlpp::drive_family<ctrlpp::estimation_demo>("estimation", ctrlpp::kHostKalmanVelocityVariance, ctrlpp::estimation_tolerance<float>());
    ctrlpp::drive_family<ctrlpp::dsp_demo>("dsp", ctrlpp::kHostBiquadOutput, ctrlpp::dsp_tolerance<float>());
    ctrlpp::drive_family<ctrlpp::trajectory_demo>("trajectory", ctrlpp::kHostTrajectoryPosition, ctrlpp::trajectory_tolerance<float>());
    ctrlpp::report_absent_predictive();
}

void run_report(demo_type &demo)
{
    ctrlpp::emit_gain(demo.K(0, 0), demo.K(0, 1));
    run_loop(demo);
    const ctrlpp::window_figures control = ctrlpp::read_window();
    run_families(control);
    const ctrlpp::golden_verdict verdict = ctrlpp::judge_golden(demo, golden_workspace);
    // ESP-IDF reports the high-water mark in bytes (vanilla FreeRTOS uses words).
    ctrlpp::emit_stack("program", static_cast<std::uint32_t>(uxTaskGetStackHighWaterMark(nullptr)), ctrlpp::kControlTaskStackBytes);
    ctrlpp::report_golden(demo, verdict, control);
    ESP_LOGI(TAG, "golden diff %s", verdict.pass ? "PASS" : "FAIL");

    // The stack probe starts only now, with the golden verdict already reported.
    // Its second pass ends in a deliberate stack overflow, which resets the chip;
    // running it concurrently would take the board down before this loop -- the
    // harness's actual witness -- had finished its steps and said whether it
    // agreed with the host.
    if constexpr(ctrlpp::kStackProbeRuns)
        ctrlpp::start_stack_probe();
}

void control_task(void *)
{
    if constexpr(!ctrlpp::kAllocSensorPresent)
        ctrlpp::stop_task_watchdog();
    report_opening();
    auto demo = demo_type::make();
    if(demo.has_value())
        run_report(*demo);
    else
        ctrlpp::report_refusal("control", ctrlpp::describe(demo.error()));
    vTaskDelete(nullptr);
}

}

extern "C" void app_main()
{
    ctrlpp::telemetry_init();
    // Eigen working set + newlib %f formatting need far more than the ~3.5 KB
    // default app_main stack.
    xTaskCreatePinnedToCore(control_task, "ctrlpp_loop", ctrlpp::kControlTaskStackBytes, nullptr, 5, nullptr, 0);
}
