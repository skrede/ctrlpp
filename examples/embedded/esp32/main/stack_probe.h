#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_H

#include <cstdint>

namespace ctrlpp {

// ESP-IDF FreeRTOS stack sizes are in bytes, not words.
inline constexpr std::uint32_t kControlTaskStackBytes = 8 * 1024;

// The probe's second pass ends in a deliberate stack overflow, which panics and
// resets the chip, so it runs in the evidence image only: the timing image then
// restarts when its capture restarts it rather than in a panic after every
// report.
inline constexpr bool kStackProbeRuns = CTRLPP_MCU_ALLOC_SENTINEL != 0;

void start_stack_probe();

}

#endif
