#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_CYCLE_COUNTER_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_CYCLE_COUNTER_H

#include <cstdint>

namespace ctrlpp {

enum class counter_state
{
    live,
    not_counting
};

const char *describe(counter_state state) noexcept;

inline constexpr std::uint32_t kCounterBits = 32;

// The other leg enables its counter and then proves it. This core's CCOUNT runs
// from reset with no trace enable and no software lock, so there is no enable
// sequence to run or to fail, and the proof is the measurement half alone. It is
// still mandatory: any state other than live means no figure may be printed.
counter_state prove_cycle_counter() noexcept;

// Each core has its own counter, so both samples of a region must come from one
// core; the timed task is pinned for that reason.
std::uint32_t cycle_sample() noexcept;

// At this part's clock the counter wraps in seconds rather than the minute the
// other leg has, so a wrap inside a timed run is a real event, not a formality.
std::uint32_t cycle_wrap_ms() noexcept;

struct cycle_span
{
    std::uint32_t cycles;
    bool wrapped;
    bool tick_crossed;
};

// An end below its start means the counter wrapped inside the region, and the
// modular difference would read as a small plausible figure.
inline cycle_span cycles_between(std::uint32_t start, std::uint32_t end) noexcept
{
    return {end - start, end < start, false};
}

// A region's opening mark carries the scheduler tick beside the cycle count.
// The tick catches what the count alone cannot: a region that ran longer than a
// whole wrap period, and a region the tick interrupt landed in. It is read
// outside the cycle samples, so the figure excludes its cost and a crossing is
// reported for a bracket a few instructions wider than the region.
struct cycle_mark
{
    std::uint32_t tick;
    std::uint32_t cycles;
};

cycle_mark open_region() noexcept;

cycle_span close_region(cycle_mark start) noexcept;

// An open and close pair around an empty region: the bias every raw figure
// carries.
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
inline void keep_value(float value) noexcept
{
    __asm__ __volatile__("" : : "f"(value));
}

}

#endif
