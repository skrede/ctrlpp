#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_DIMENSION_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_STACK_PROBE_DIMENSION_H

#include "freertos/FreeRTOS.h"
#include "freertos/semphr.h"

#include <cstddef>
#include <cstdint>

namespace ctrlpp::probe
{

inline constexpr const char* probe_tag = "ctrlpp_esp32";

// Returns once the dimension's own task, given stack bytes, has reported. Each dimension below is
// instantiated in its own translation unit, generated from stack_probe_dimension.cpp.in by
// main/CMakeLists.txt, since one unit holding every dimension needed 11 GB to compile for Xtensa; a
// dimension declared here but not generated there fails the link.
template <std::size_t NX, std::size_t NU>
void run_one(std::uint32_t stack, SemaphoreHandle_t done);

extern template void run_one<2, 1>(std::uint32_t, SemaphoreHandle_t);
extern template void run_one<3, 1>(std::uint32_t, SemaphoreHandle_t);
extern template void run_one<4, 2>(std::uint32_t, SemaphoreHandle_t);
extern template void run_one<5, 1>(std::uint32_t, SemaphoreHandle_t);
extern template void run_one<6, 2>(std::uint32_t, SemaphoreHandle_t);
extern template void run_one<8, 2>(std::uint32_t, SemaphoreHandle_t);

}

#endif
