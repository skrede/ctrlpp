#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TELEMETRY_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TELEMETRY_H

#include "alloc_sensor.h"
#include "cycle_record.h"
#include "cycle_counter.h"

#include "telemetry_format.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

inline constexpr std::size_t kLineBytes = 320;

void telemetry_init();

// Every line opens with exactly one tag, so a block is the set of lines carrying
// its tag and a missing class reads as a missing tag.
void send_line(stream_tag tag, const char *body);

void emit_identity();

void emit_allocation_canary(const canary_verdict &verdict);

void emit_counter_canary(counter_state state);

void emit_header(claim_class claim);

inline constexpr std::size_t kExactValues = column_count(claim_class::exact) - 1;

void emit_exact_row(std::int32_t k, const std::array<float, kExactValues> &values);

// The figures are raw: each carries the instrument overhead, which the report
// publishes on its own rows. A region whose every sample wrapped prints a zero
// count and no figures.
void emit_spread_row(const timing_region &region, const cycle_distribution &d);

// How many of the region's samples the scheduler tick interrupt landed in: a
// lower bound on how often it was interrupted, since other interrupts leave no
// mark the instrument can read.
void emit_interruption(const timing_region &region, const cycle_distribution &d);

void emit_gain(float k0, float k1);

void emit_stack(const char *region, std::uint32_t free_bytes, std::uint32_t stack_bytes);

}

#endif
