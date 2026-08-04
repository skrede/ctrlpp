/// ctrlpp ESP32 on-device control-loop leg. Drives the shared float control
/// kernel from a fixed-cadence FreeRTOS task on real Xtensa silicon: designs the
/// infinite-horizon LQR gain once, runs the closed loop at 50 Hz, streams the
/// trajectory as CSV over UART2 (telemetry adapter, ESP TX GPIO17 -> USB RX), and
/// diffs the on-target float result against an independent host-double golden.

#include "golden_reference.h"
#include "control_loop_demo.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "freertos/semphr.h"

#include "esp_log.h"
#include "driver/uart.h"
#include "esp_rom_sys.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cstddef>

namespace
{

constexpr const char* TAG = "ctrlpp_esp32";

// Telemetry UART (UART2): ESP TX GPIO17 -> USB RX, ESP RX GPIO16 <- USB TX. GPIO
// 16/17 are free on WROOM but the PSRAM lines on WROVER -- keep them named.
constexpr uart_port_t kTelemetryUart = UART_NUM_2;
constexpr int         kTxPin         = 17;
constexpr int         kRxPin         = 16;
constexpr int         kBaud          = 115200;

void telemetry_init()
{
    uart_config_t cfg{};
    cfg.baud_rate  = kBaud;
    cfg.data_bits  = UART_DATA_8_BITS;
    cfg.parity     = UART_PARITY_DISABLE;
    cfg.stop_bits  = UART_STOP_BITS_1;
    cfg.flow_ctrl  = UART_HW_FLOWCTRL_DISABLE;
    cfg.source_clk = UART_SCLK_DEFAULT;

    uart_driver_install(kTelemetryUart, 512, 0, 0, nullptr, 0);
    uart_param_config(kTelemetryUart, &cfg);
    uart_set_pin(kTelemetryUart, kTxPin, kRxPin, UART_PIN_NO_CHANGE, UART_PIN_NO_CHANGE);
}

void telemetry_write(const char* line)
{
    uart_write_bytes(kTelemetryUart, line, std::strlen(line));
}

void start_stack_probe();

void control_task(void*)
{
    auto demo = ctrlpp::control_loop_demo<float>::make();
    if(!demo.has_value())
    {
        ESP_LOGE(TAG, "lqr_gain REFUSED the plant on device -- %s", ctrlpp::describe(demo.error()));
        return;
    }

    const double dev0 = static_cast<double>(demo->K(0, 0)) - ctrlpp::kHostK0;
    const double dev1 = static_cast<double>(demo->K(0, 1)) - ctrlpp::kHostK1;
    ESP_LOGI(TAG, "gain    K = [%.7f, %.7f]", static_cast<double>(demo->K(0, 0)),
             static_cast<double>(demo->K(0, 1)));
    ESP_LOGI(TAG, "host    K = [%.7f, %.7f]", ctrlpp::kHostK0, ctrlpp::kHostK1);
    ESP_LOGI(TAG, "float-vs-host gain dev = [%.3e, %.3e]", dev0, dev1);

    telemetry_write("k,t_s,x0,x1,u\n");

    TickType_t       next   = xTaskGetTickCount();
    const TickType_t period = pdMS_TO_TICKS(20);   // 50 Hz cadence
    char buf[80];
    for(int k = 0; k <= ctrlpp::kSteps; ++k)
    {
        const float u = demo->step();
        std::snprintf(buf, sizeof(buf), "%d,%.4f,%.6f,%.6f,%.6f\n", k,
                      static_cast<double>(k * ctrlpp::kDt),
                      static_cast<double>(demo->x[0]),
                      static_cast<double>(demo->x[1]), static_cast<double>(u));
        telemetry_write(buf);
        vTaskDelayUntil(&next, period);
    }

    // ESP-IDF reports the high-water mark in bytes (vanilla FreeRTOS uses words).
    ESP_LOGI(TAG, "stack high-water = %u bytes",
             static_cast<unsigned>(uxTaskGetStackHighWaterMark(nullptr)));

    const double final_norm = static_cast<double>(demo->x.norm());
    const double err        = std::fabs(final_norm - ctrlpp::kHostFinalNorm);
    ESP_LOGI(TAG, "settling: final |x| = %.3e (host = %.3e, err = %.3e)",
             final_norm, ctrlpp::kHostFinalNorm, err);
    ESP_LOGI(TAG, "golden diff %s (tol = %.1e)",
             err < ctrlpp::kEsp32FloatTol ? "PASS" : "FAIL", ctrlpp::kEsp32FloatTol);

    // The stack probe starts only now, with the golden verdict already reported.
    // Its second pass ends in a deliberate stack overflow, which resets the chip;
    // running it concurrently would take the board down before this loop -- the
    // harness's actual witness -- had finished its steps and said whether it
    // agreed with the host.
    start_stack_probe();

    vTaskDelete(nullptr);
}

// ---------------------------------------------------------------------------
// Stack probe: how large a state dimension the discrete solve fits in this task
// ---------------------------------------------------------------------------
//
// Every stack figure this project publishes for the Riccati path so far comes
// from a host compiler's -fstack-usage report. That is a PREDICTION about a
// different architecture, a different calling convention and a different
// register file. This probe is the test.
//
// Two passes, because "what does this dimension cost" and "does this dimension
// fit" are different questions and one sweep cannot answer both.
//
// PASS 1 gives every dimension a generous stack, so each one returns and
// reports what it actually needs. PASS 2 gives each dimension the control
// task's own stack size, so the figures share that denominator and the first
// dimension that does not fit establishes the supported maximum.
//
// Each dimension runs in ITS OWN task. The high-water mark is a running
// minimum over a task's life, so one task sweeping every dimension would report
// the deepest dimension's figure against all of them. One task per dimension
// makes each figure that dimension's own.
//
// Dimensions run in INCREASING order and each reports before the next is
// attempted, so a dimension that does not fit leaves a complete record of every
// dimension below it.

constexpr std::uint32_t kTaskStackBytes = 8 * 1024;

/// Deliberately generous, so that pass one measures what each dimension NEEDS
/// rather than only whether it fitted the control task's stack. A dimension that
/// overflows here would be reported the same way, by name.
constexpr std::uint32_t kHeadroomStackBytes = 48 * 1024;

/// The discrete damped chain, matching the host benchmark's corpus so the two
/// describe the same family: a forward-Euler step of a chain with -0.5 on the
/// diagonal and 1.0 on the superdiagonal. Q = I and R = 0.1 I put the weight
/// scale at exactly one, so the entry point runs its acceptance check once.
template <std::size_t NX, std::size_t NU>
void probe_dimension(std::uint32_t stack)
{
    constexpr int   n  = int(NX);
    constexpr int   nu = int(NU);
    constexpr float dt = 0.01f;

    Eigen::Matrix<float, n, n> A = Eigen::Matrix<float, n, n>::Identity();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5f;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt;

    Eigen::Matrix<float, n, nu> B     = Eigen::Matrix<float, n, nu>::Zero();
    const std::size_t           group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
        B(int((j + 1) * group - 1), int(j)) = dt;

    Eigen::Matrix<float, n, n>   Q = Eigen::Matrix<float, n, n>::Identity();
    Eigen::Matrix<float, nu, nu> R = 0.1f * Eigen::Matrix<float, nu, nu>::Identity();

    auto solved = ctrlpp::dare<float, NX, NU>(A, B, Q, R);

    // Reported on the CONSOLE endpoint. A task that has overflowed cannot be
    // trusted to complete a telemetry write, and this line is the evidence.
    //
    // The high-water mark is the MINIMUM free this task ever had, so reading it
    // from a task that solved exactly one dimension makes the figure that
    // dimension's own requirement. A single task sweeping every dimension would
    // report a running minimum instead, which answers "did the sweep fit" and
    // not "what does this dimension cost".
    const unsigned free_bytes = static_cast<unsigned>(uxTaskGetStackHighWaterMark(nullptr));
    ESP_LOGI(TAG, "probe NX=%u NU=%u: %s, used %u of %u bytes (high-water %u free)",
             static_cast<unsigned>(NX), static_cast<unsigned>(NU),
             solved.has_value() ? "solved" : ctrlpp::describe(solved.error()),
             static_cast<unsigned>(stack) - free_bytes, static_cast<unsigned>(stack),
             free_bytes);
}

/// One task per dimension, so each figure is that dimension's own requirement.
template <std::size_t NX, std::size_t NU>
void run_one(std::uint32_t stack, SemaphoreHandle_t done)
{
    struct arg_t { std::uint32_t stack; SemaphoreHandle_t done; };
    static arg_t arg;
    arg = {stack, done};
    auto body = [](void* raw) {
        auto* a = static_cast<arg_t*>(raw);
        probe_dimension<NX, NU>(a->stack);
        xSemaphoreGive(a->done);
        vTaskDelete(nullptr);
    };
    char name[16];
    std::snprintf(name, sizeof(name), "probe_nx%u", static_cast<unsigned>(NX));
    if(xTaskCreatePinnedToCore(body, name, stack, &arg, 4, nullptr, 0) != pdPASS)
    {
        ESP_LOGE(TAG, "probe NX=%u: could not create a task with a %u byte stack",
                 static_cast<unsigned>(NX), static_cast<unsigned>(stack));
        return;
    }
    xSemaphoreTake(done, portMAX_DELAY);
}

/// Pass one: a stack large enough that every dimension survives, so each one
/// reports what it actually needs rather than only whether it fitted.
void headroom_pass(SemaphoreHandle_t done)
{
    ESP_LOGI(TAG, "PASS 1 -- requirement per dimension, %u byte stack, scalar = float",
             static_cast<unsigned>(kHeadroomStackBytes));
    run_one<2, 1>(kHeadroomStackBytes, done);
    run_one<3, 1>(kHeadroomStackBytes, done);
    run_one<4, 2>(kHeadroomStackBytes, done);
    run_one<5, 1>(kHeadroomStackBytes, done);
    run_one<6, 2>(kHeadroomStackBytes, done);
    run_one<8, 2>(kHeadroomStackBytes, done);
    ESP_LOGI(TAG, "PASS 1 complete -- every dimension returned");
}

/// Pass two: the control task's own stack size. This one is expected to end in a
/// named overflow, and that is the point -- it establishes the supported maximum
/// by reaching the first dimension that does not fit.
void limit_pass(SemaphoreHandle_t done)
{
    ESP_LOGI(TAG, "PASS 2 -- supported maximum at the control task's own %u byte stack",
             static_cast<unsigned>(kTaskStackBytes));
    run_one<2, 1>(kTaskStackBytes, done);
    run_one<3, 1>(kTaskStackBytes, done);
    run_one<4, 2>(kTaskStackBytes, done);
    run_one<5, 1>(kTaskStackBytes, done);
    run_one<6, 2>(kTaskStackBytes, done);
    run_one<8, 2>(kTaskStackBytes, done);
    ESP_LOGI(TAG, "PASS 2 complete -- no dimension overflowed, which was not expected");
}

void stack_probe_task(void*)
{
    SemaphoreHandle_t done = xSemaphoreCreateBinary();
    headroom_pass(done);
    limit_pass(done);
    vTaskDelete(nullptr);
}

/// The driver task itself needs almost no stack: every solve happens in a task
/// this one creates, and it only waits on them.
void start_stack_probe()
{
    xTaskCreatePinnedToCore(stack_probe_task, "ctrlpp_probe", 3 * 1024, nullptr, 4, nullptr, 0);
}

}

/// Name the task that overflowed instead of letting the board reset unexplained.
/// The framework supplies a weak default; this overrides it so a negative result
/// is evidence rather than an absence of output.
///
/// Both a log line and a ROM print: the log path takes locks and allocates, and
/// a task that has just overrun its stack is the last context that should be
/// trusted to complete either. `esp_rom_printf` writes the console UART
/// directly and is what survives if the log path does not.
extern "C" void vApplicationStackOverflowHook(TaskHandle_t /*task*/, char* name)
{
    esp_rom_printf("STACK OVERFLOW in task '%s' -- the dimension above the last "
                   "reported probe line does not fit %u bytes\n",
                   name != nullptr ? name : "<unnamed>",
                   static_cast<unsigned>(kTaskStackBytes));
}

extern "C" void app_main()
{
    telemetry_init();
    // Eigen working set + newlib %f formatting need far more than the ~3.5 KB
    // default app_main stack; ESP-IDF FreeRTOS stack sizes are in BYTES.
    xTaskCreatePinnedToCore(control_task, "ctrlpp_loop", kTaskStackBytes, nullptr, 5, nullptr, 0);
    // The probe is started by control_task once its golden comparison has been
    // reported, NOT here. Pass two ends in a deliberate stack overflow, which
    // resets the chip; starting the probe concurrently would take the board down
    // before the control loop -- this harness's actual witness -- had finished
    // its 201 steps and printed its verdict.
}
