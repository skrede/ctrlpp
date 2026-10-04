#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_CYCLE_COUNTER_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_CYCLE_COUNTER_H

#include <cstdint>

namespace ctrlpp {

enum class counter_state
{
    live,
    absent,
    not_enabled,
    not_counting
};

const char *describe(counter_state state) noexcept;

// Runs the enable sequence from a forced-off state and then proves the counter
// is running; any state other than live means no figure may be printed.
counter_state enable_cycle_counter() noexcept;

// The measurement half of the proof alone, for a caller that needs to know the
// counter is still running without enabling it again.
counter_state prove_cycle_counter() noexcept;

// Out of line, so no caller has to parse the device header: the call costs a
// few cycles, and measure_overhead charges them to the published overhead.
std::uint32_t cycle_sample() noexcept;

struct cycle_span
{
    std::uint32_t cycles;
    bool wrapped;
};

// An end below its start means the counter wrapped inside the region, and the
// modular difference would read as a small plausible figure.
inline cycle_span cycles_between(std::uint32_t start, std::uint32_t end) noexcept
{
    return {end - start, end < start};
}

// A sample pair around an empty region: the bias every raw figure carries.
cycle_span measure_overhead() noexcept;

// Forces an object into memory with its stores complete, so the work that
// produced it cannot be removed or moved past the closing sample.
template<class T>
void keep(const T &object) noexcept
{
    __asm__ __volatile__("" : : "r"(&object) : "memory");
}

// The same for a value held in a floating-point register, at no cost beyond
// computing it.
inline void keep_value(double value) noexcept
{
    __asm__ __volatile__("" : : "w"(value));
}

}

#endif
