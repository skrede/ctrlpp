// ctrlpp ESP32 on-device leg (FreeRTOS, single-precision FPU).
//
// Drives the shared float control kernel from a fixed-cadence FreeRTOS task on
// real Xtensa silicon: prints the reset that started the run and the posture it
// is taken under, designs the infinite-horizon LQR gain once, runs the closed
// loop at 50 Hz while streaming the trajectory's exact rows over UART2
// (telemetry adapter, ESP TX GPIO5 -> USB RX), runs the estimation,
// signal-processing and trajectory families, and diffs every family's float
// result against its host-double golden.

#include "telemetry.h"
#include "stack_probe.h"
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

void report_opening()
{
    ctrlpp::emit_identity();
    ctrlpp::emit_posture(ctrlpp::read_reset_cause());
    ctrlpp::emit_header(ctrlpp::claim_class::spread);
    ctrlpp::emit_canaries();
}

void run_loop(demo_type &demo)
{
    ctrlpp::emit_header(ctrlpp::claim_class::exact);
    TickType_t next         = xTaskGetTickCount();
    const TickType_t period = pdMS_TO_TICKS(20);
    for(std::int32_t k = 0; k <= ctrlpp::kSteps; ++k)
    {
        const float u = demo.step();
        ctrlpp::emit_exact_row(k, {static_cast<float>(k * ctrlpp::kDt), demo.x[0], demo.x[1], u});
        vTaskDelayUntil(&next, period);
    }
}

void run_families()
{
    ctrlpp::report_unmeasured_allocations();
    ctrlpp::drive_family<ctrlpp::estimation_demo>("estimation", ctrlpp::kHostKalmanVelocityVariance, ctrlpp::estimation_tolerance<float>());
    ctrlpp::drive_family<ctrlpp::dsp_demo>("dsp", ctrlpp::kHostBiquadOutput, ctrlpp::dsp_tolerance<float>());
    ctrlpp::drive_family<ctrlpp::trajectory_demo>("trajectory", ctrlpp::kHostTrajectoryPosition, ctrlpp::trajectory_tolerance<float>());
    ctrlpp::report_absent_predictive();
}

void run_report(demo_type &demo)
{
    ctrlpp::emit_gain(demo.K(0, 0), demo.K(0, 1));
    run_loop(demo);
    run_families();
    const ctrlpp::golden_verdict verdict = ctrlpp::judge_golden(demo, golden_workspace);
    // ESP-IDF reports the high-water mark in bytes (vanilla FreeRTOS uses words).
    ctrlpp::emit_stack("program", static_cast<std::uint32_t>(uxTaskGetStackHighWaterMark(nullptr)), ctrlpp::kControlTaskStackBytes);
    ctrlpp::report_golden(demo, verdict);
    ESP_LOGI(TAG, "golden diff %s", verdict.pass ? "PASS" : "FAIL");

    // The stack probe starts only now, with the golden verdict already reported.
    // Its second pass ends in a deliberate stack overflow, which resets the chip;
    // running it concurrently would take the board down before this loop -- the
    // harness's actual witness -- had finished its steps and said whether it
    // agreed with the host.
    ctrlpp::start_stack_probe();
}

void control_task(void *)
{
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
