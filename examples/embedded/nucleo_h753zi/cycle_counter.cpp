#include "cycle_counter.h"

#include "stm32h7xx.h"

#include <cstdint>

namespace {

// The CoreSight software-lock key (Arm CoreSight Architecture Specification,
// Lock Access Register). The core-support headers carry the register but not
// the key.
constexpr std::uint32_t kCoreSightUnlockKey = 0xC5ACCE55U;
constexpr std::uint32_t kCoreSightLock      = 0U;

// The trace block and the counter sit in the debug power domain, which a system
// reset does not clear, so a previous image or a debugger can leave the counter
// running. Turning it off first makes the proof depend on this image's sequence
// alone.
void force_counter_off() noexcept
{
    DWT->LAR = kCoreSightUnlockKey;
    DWT->CTRL &= ~DWT_CTRL_CYCCNTENA_Msk;
    DWT->LAR = kCoreSightLock;
    DCB->DEMCR &= ~DCB_DEMCR_TRCENA_Msk;
    __DSB();
}

// The unlock is written unconditionally rather than gated on the lock-status
// register: the core-support headers define no masks for that register, and
// the Arm Cortex-M7 Software Developers Errata Notice records that its status
// bit reads as one for debugger accesses. Liveness is proved by measurement.
void run_enable_sequence() noexcept
{
    DCB->DEMCR |= DCB_DEMCR_TRCENA_Msk;
    DWT->LAR    = kCoreSightUnlockKey;
    DWT->CYCCNT = 0U;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
    __DSB();
}

}

namespace ctrlpp {

const char *describe(counter_state state) noexcept
{
    switch(state)
    {
        case counter_state::live:
            return "live";
        case counter_state::absent:
            return "the core reports no cycle counter";
        case counter_state::not_enabled:
            return "the counter enable did not read back set";
        case counter_state::not_counting:
            return "two samples around a few no-operations did not differ";
    }
    return "unknown";
}

// The barriers keep the read from being moved across the region it brackets:
// the M7 issues in order but completes out of order, and to the compiler the
// register read is an ordinary load apart from its volatility.
std::uint32_t cycle_sample() noexcept
{
    __DSB();
    __ISB();
    const std::uint32_t value = DWT->CYCCNT;
    __DSB();
    __ISB();
    return value;
}

counter_state enable_cycle_counter() noexcept
{
    force_counter_off();
    run_enable_sequence();
    return prove_cycle_counter();
}

counter_state prove_cycle_counter() noexcept
{
    const std::uint32_t control = DWT->CTRL;
    if((control & DWT_CTRL_NOCYCCNT_Msk) != 0U)
        return counter_state::absent;
    if((control & DWT_CTRL_CYCCNTENA_Msk) == 0U)
        return counter_state::not_enabled;

    const std::uint32_t before = cycle_sample();
    __NOP();
    __NOP();
    __NOP();
    __NOP();
    const cycle_span span = cycles_between(before, cycle_sample());
    return !span.wrapped && span.cycles > 0U ? counter_state::live : counter_state::not_counting;
}

cycle_span measure_overhead() noexcept
{
    const std::uint32_t start = cycle_sample();
    return cycles_between(start, cycle_sample());
}

}
