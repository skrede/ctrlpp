#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_H

#include <cstdint>

namespace ctrlpp {

// ESP-IDF FreeRTOS stack sizes are in bytes, not words.
inline constexpr std::uint32_t kControlTaskStackBytes = 8 * 1024;

// Its second pass ends in a deliberate stack overflow, which resets the chip.
void start_stack_probe();

}

#endif
