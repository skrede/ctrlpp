#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_DIMENSION_IMPL_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_DIMENSION_IMPL_H

#include "stack_probe_dimension.h"

#include "control_loop_demo.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "freertos/semphr.h"

#include "esp_log.h"

#include <Eigen/Dense>

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp::probe {

// The discrete damped chain, matching the host benchmark's corpus so the two
// describe the same family: a forward-Euler step of a chain with -0.5 on the
// diagonal and 1.0 on the superdiagonal. Q = I and R = 0.1 I put the weight
// scale at exactly one, so the entry point runs its acceptance check once.
template<std::size_t NX, std::size_t NU>
const char *solve_chain()
{
    constexpr int n    = int(NX);
    constexpr int nu   = int(NU);
    constexpr float dt = 0.01f;

    Eigen::Matrix<float, n, n> A = Eigen::Matrix<float, n, n>::Identity();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5f;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt;

    Eigen::Matrix<float, n, nu> B = Eigen::Matrix<float, n, nu>::Zero();
    const std::size_t group       = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
        B(int((j + 1) * group - 1), int(j)) = dt;

    Eigen::Matrix<float, n, n> Q   = Eigen::Matrix<float, n, n>::Identity();
    Eigen::Matrix<float, nu, nu> R = 0.1f * Eigen::Matrix<float, nu, nu>::Identity();

    const auto solved = ctrlpp::dare<float, NX, NU>(A, B, Q, R);
    return solved.has_value() ? "solved" : ctrlpp::describe(solved.error());
}

// Reported on the CONSOLE endpoint. A task that has overflowed cannot be
// trusted to complete a telemetry write, and this line is the evidence.
// The high-water mark is the MINIMUM free this task ever had, so reading it
// from a task that solved exactly one dimension makes the figure that
// dimension's own requirement. A single task sweeping every dimension would
// report a running minimum instead, which answers "did the sweep fit" and
// not "what does this dimension cost".
template<std::size_t NX, std::size_t NU>
void probe_dimension(std::uint32_t stack)
{
    const char *outcome            = solve_chain<NX, NU>();
    const std::uint32_t free_bytes = static_cast<std::uint32_t>(uxTaskGetStackHighWaterMark(nullptr));
    ESP_LOGI(kProbeTag, "probe NX=%" PRIu32 " NU=%" PRIu32 ": %s, used %" PRIu32 " of %" PRIu32 " bytes (high-water %" PRIu32 " free)", static_cast<std::uint32_t>(NX),
             static_cast<std::uint32_t>(NU), outcome, stack - free_bytes, stack, free_bytes);
}

struct probe_arg
{
    std::uint32_t stack;
    SemaphoreHandle_t done;
};

template<std::size_t NX, std::size_t NU>
void run_one(std::uint32_t stack, SemaphoreHandle_t done)
{
    static probe_arg arg;
    arg       = {stack, done};
    auto body = [](void *raw)
    {
        auto *a = static_cast<probe_arg *>(raw);
        probe_dimension<NX, NU>(a->stack);
        xSemaphoreGive(a->done);
        vTaskDelete(nullptr);
    };
    char name[16];
    std::snprintf(name, sizeof(name), "probe_nx%" PRIu32, static_cast<std::uint32_t>(NX));
    if(xTaskCreatePinnedToCore(body, name, stack, &arg, 4, nullptr, 0) != pdPASS)
    {
        ESP_LOGE(kProbeTag, "probe NX=%" PRIu32 ": could not create a task with a %" PRIu32 " byte stack", static_cast<std::uint32_t>(NX), stack);
        return;
    }
    xSemaphoreTake(done, portMAX_DELAY);
}

}

#endif
