#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_TELEMETRY_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_TELEMETRY_H

#include "alloc_sensor.h"
#include "cycle_record.h"
#include "cycle_counter.h"

#include "telemetry_format.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

// Every line of the report opens with exactly one tag, so a block is the set of
// lines carrying its tag and a missing class reads as a missing tag.
void open_line(stream_tag tag);

void emit_identity();

void emit_allocation_canary(const canary_verdict &verdict);

void emit_counter_canary(counter_state state);

// Lists the class's fields in their table order; every row of the class that
// follows carries its values in the same order.
void emit_header(claim_class claim);

inline constexpr std::size_t kExactValues = column_count(claim_class::exact) - 1;

// One row per loop index, in index order, so two runs compare line by line and a
// reordering is a difference rather than a wash.
void emit_exact_row(std::int32_t k, const std::array<double, kExactValues> &values);

// The figures are raw: each carries the instrument overhead, which the report
// publishes on its own rows. A region whose every sample wrapped prints a zero
// count and no figures.
void emit_spread_row(const timing_region &region, const cycle_distribution &d);

void emit_setup_heap();

void emit_gain(double k0, double k1);

}

#endif
