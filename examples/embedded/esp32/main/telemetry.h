#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TELEMETRY_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TELEMETRY_H

#include "reset_reason.h"

#include "telemetry_format.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

inline constexpr std::size_t kLineBytes = 256;

void telemetry_init();

// Every line opens with exactly one tag, so a block is the set of lines carrying
// its tag and a missing class reads as a missing tag.
void send_line(stream_tag tag, const char *body);

void emit_identity();

void emit_posture(const reset_cause &cause);

void emit_canaries();

void emit_header(claim_class claim);

inline constexpr std::size_t kExactValues = column_count(claim_class::exact) - 1;

void emit_exact_row(std::int32_t k, const std::array<float, kExactValues> &values);

void emit_gain(float k0, float k1);

void emit_stack(const char *region, std::uint32_t free_bytes, std::uint32_t stack_bytes);

}

#endif
