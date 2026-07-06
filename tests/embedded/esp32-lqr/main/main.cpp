/// ctrlpp ESP32 on-device LQR check. Designs an infinite-horizon LQR gain for a
/// discrete double-integrator (50 Hz) with `float` scalars on real silicon,
/// runs the closed loop from a unit position error, and:
///   - streams the trajectory as CSV over UART2 (the telemetry adapter,
///     ESP TX GPIO17 -> USB RX) so it can be captured and diffed on the host;
///   - logs a human-readable summary over the console UART (the flash port),
///     including the on-device float gain, its deviation from the host-double
///     reference, closed-loop settling, and per-step cycle/microsecond timing.
///
/// The reference constants below are produced by a host-double build of the
/// exact same problem, following the "never trust the self-report -- diff
/// against an independent reference" discipline. Exercises
/// ctrlpp::lqr_gain (DARE-backed) with exceptions and RTTI off.

#include "ctrlpp/control.h"

#include "driver/uart.h"
#include "esp_timer.h"
#include "esp_cpu.h"
#include "esp_log.h"

#include <Eigen/Dense>

#include <cstdio>
#include <cstring>
#include <cstdint>

namespace
{

constexpr const char* TAG = "ctrlpp_lqr";

using scalar = float;

// --- Telemetry UART (UART2): ESP TX GPIO17 -> USB RX, ESP RX GPIO16 <- USB TX.
constexpr uart_port_t kTelemetryUart = UART_NUM_2;
constexpr int         kTxPin         = 17;
constexpr int         kRxPin         = 16;
constexpr int         kBaud          = 115200;

// --- Discrete double-integrator at 50 Hz.
constexpr scalar kDt    = 0.02f;
constexpr int    kSteps = 200;

// --- Host-double reference for this exact problem (independent cross-check).
constexpr double kHostK0        = 9.467093672;
constexpr double kHostK1        = 5.281739638;
constexpr double kHostFinalNorm = 2.020686726e-05;

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

}

extern "C" void app_main()
{
    telemetry_init();

    // --- Plant: discrete double-integrator, u accelerates the second state. ---
    Eigen::Matrix<scalar, 2, 2> A;
    A << scalar{1}, kDt, scalar{0}, scalar{1};
    Eigen::Matrix<scalar, 2, 1> B;
    B << scalar{0.5} * kDt * kDt, kDt;
    Eigen::Matrix<scalar, 2, 2> Q;
    Q << scalar{10}, scalar{0}, scalar{0}, scalar{1};
    Eigen::Matrix<scalar, 1, 1> R;
    R << scalar{0.1};

    // --- Design the infinite-horizon gain on-device (float DARE). ---
    const auto gain = ctrlpp::lqr_gain<scalar, 2, 1>(A, B, Q, R);
    if (!gain.has_value())
    {
        ESP_LOGE(TAG, "lqr_gain FAILED on device -- Riccati did not converge");
        return;
    }
    const Eigen::Matrix<scalar, 1, 2> K = *gain;

    const double dev0 = static_cast<double>(K(0, 0)) - kHostK0;
    const double dev1 = static_cast<double>(K(0, 1)) - kHostK1;
    ESP_LOGI(TAG, "gain    K = [%.7f, %.7f]", static_cast<double>(K(0, 0)),
             static_cast<double>(K(0, 1)));
    ESP_LOGI(TAG, "host    K = [%.7f, %.7f]", kHostK0, kHostK1);
    ESP_LOGI(TAG, "float-vs-host gain dev = [%.3e, %.3e]", dev0, dev1);

    // --- Per-step feedback-law timing (K*x), the RT-critical inner product. ---
    Eigen::Vector<scalar, 2> probe;
    probe << scalar{0.3}, scalar{-0.1};
    volatile scalar sink = 0.f;
    for (int i = 0; i < 16; ++i)  // warm caches
    {
        sink += (K * probe)(0);
    }
    constexpr int law_iters = 2000;
    const uint32_t c0 = esp_cpu_get_cycle_count();
    const int64_t  t0 = esp_timer_get_time();
    for (int i = 0; i < law_iters; ++i)
    {
        probe[0] = sink;  // break common-subexpression elimination
        sink += (K * probe)(0);
    }
    const int64_t  t1 = esp_timer_get_time();
    const uint32_t c1 = esp_cpu_get_cycle_count();
    ESP_LOGI(TAG, "feedback law: %.4f us/eval, %lu cycles/eval (avg of %d)",
             static_cast<double>(t1 - t0) / law_iters,
             static_cast<unsigned long>((c1 - c0) / law_iters), law_iters);

    // --- Closed loop from a unit position error; stream CSV over UART2. ---
    telemetry_write("k,t_s,x0,x1,u\n");
    Eigen::Vector<scalar, 2> x;
    x << scalar{1}, scalar{0};
    scalar u_peak = 0.f;
    char buf[80];
    for (int k = 0; k <= kSteps; ++k)
    {
        const scalar u = -(K * x)(0);
        if (std::fabs(u) > std::fabs(u_peak))
        {
            u_peak = u;
        }
        std::snprintf(buf, sizeof(buf), "%d,%.4f,%.6f,%.6f,%.6f\n", k,
                      static_cast<double>(k * kDt), static_cast<double>(x[0]),
                      static_cast<double>(x[1]), static_cast<double>(u));
        telemetry_write(buf);
        x = A * x + B * u;
    }

    const double final_norm = static_cast<double>(x.norm());
    ESP_LOGI(TAG, "closed loop: %d steps @ %.0f Hz, peak |u| = %.4f", kSteps,
             1.0 / static_cast<double>(kDt), static_cast<double>(u_peak));
    ESP_LOGI(TAG, "settling: final |x| = %.3e (host = %.3e, ratio = %.3f)",
             final_norm, kHostFinalNorm, final_norm / kHostFinalNorm);
    ESP_LOGI(TAG, "done -- CSV trajectory streamed on UART2 (GPIO%d)", kTxPin);

    (void)sink;
}
