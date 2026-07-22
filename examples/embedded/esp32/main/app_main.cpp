/// ctrlpp ESP32 on-device control-loop leg. Drives the shared float control
/// kernel from a fixed-cadence FreeRTOS task on real Xtensa silicon: designs the
/// infinite-horizon LQR gain once, runs the closed loop at 50 Hz, streams the
/// trajectory as CSV over UART2 (telemetry adapter, ESP TX GPIO17 -> USB RX), and
/// diffs the on-target float result against an independent host-double golden.

#include "golden_reference.h"
#include "control_loop_demo.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_log.h"
#include "driver/uart.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>
#include <cstring>

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

void control_task(void*)
{
    auto demo = ctrlpp::control_loop_demo<float>::make();
    if(!demo.has_value())
    {
        ESP_LOGE(TAG, "lqr_gain FAILED on device -- Riccati did not converge");
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

    vTaskDelete(nullptr);
}

}

extern "C" void app_main()
{
    telemetry_init();
    // Eigen working set + newlib %f formatting need far more than the ~3.5 KB
    // default app_main stack; ESP-IDF FreeRTOS stack sizes are in BYTES.
    xTaskCreatePinnedToCore(control_task, "ctrlpp_loop", 8 * 1024, nullptr, 5, nullptr, 0);
}
