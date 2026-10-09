#include "cycle_counter.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_cpu.h"
#include "esp_rom_sys.h"

#include <cstdint>
#include <algorithm>

namespace {

constexpr std::uint64_t kCounterStates         = std::uint64_t{1} << ctrlpp::kCounterBits;
constexpr std::uint64_t kMicrosecondsPerMs     = 1000;
constexpr std::uint64_t kMillisecondsPerSecond = 1000;

static_assert(sizeof(TickType_t) == sizeof(std::uint32_t), "a mark stores the tick count in 32 bits");

std::uint32_t tick_now() noexcept
{
    return static_cast<std::uint32_t>(xTaskGetTickCount());
}

// A region at least one wrap period long crosses at least this many ticks, so
// the comparison cannot miss a whole-period wrap; it can only flag a region a
// tick short of one.
std::uint32_t wrap_ticks() noexcept
{
    return static_cast<std::uint32_t>(std::uint64_t{ctrlpp::cycle_wrap_ms()} * configTICK_RATE_HZ / kMillisecondsPerSecond);
}

}

namespace ctrlpp {

const char *describe(counter_state state) noexcept
{
    switch(state)
    {
        case counter_state::live:
            return "live";
        case counter_state::not_counting:
            return "two samples around a few no-operations did not differ";
    }
    return "unknown";
}

// The core executes in order, so the reordering risk is the compiler's; the
// memory clobbers keep loads and stores from moving across the read.
std::uint32_t cycle_sample() noexcept
{
    __asm__ __volatile__("" : : : "memory");
    const std::uint32_t value = esp_cpu_get_cycle_count();
    __asm__ __volatile__("" : : : "memory");
    return value;
}

std::uint32_t cycle_wrap_ms() noexcept
{
    const std::uint64_t ticks_per_us = std::max<std::uint32_t>(esp_rom_get_cpu_ticks_per_us(), 1U);
    return static_cast<std::uint32_t>(kCounterStates / (ticks_per_us * kMicrosecondsPerMs));
}

counter_state prove_cycle_counter() noexcept
{
    const std::uint32_t before = cycle_sample();
    __asm__ __volatile__("nop\n\tnop\n\tnop\n\tnop");
    const cycle_span span = cycles_between(before, cycle_sample());
    return !span.wrapped && span.cycles > 0U ? counter_state::live : counter_state::not_counting;
}

cycle_mark open_region() noexcept
{
    const std::uint32_t tick = tick_now();
    return {tick, cycle_sample()};
}

cycle_span close_region(cycle_mark start) noexcept
{
    const std::uint32_t end     = cycle_sample();
    const std::uint32_t elapsed = tick_now() - start.tick;
    cycle_span span             = cycles_between(start.cycles, end);
    span.wrapped                = span.wrapped || elapsed >= wrap_ticks();
    span.tick_crossed           = elapsed > 0U;
    return span;
}

cycle_span measure_overhead() noexcept
{
    return close_region(open_region());
}

}
