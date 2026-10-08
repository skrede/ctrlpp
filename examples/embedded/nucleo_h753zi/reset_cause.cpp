#include "reset_cause.h"

#include "stm32h7xx.h"

#include <array>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp {

namespace {

struct named_flag
{
    std::uint32_t mask;
    const char *name;
};

struct event_signature
{
    std::uint32_t flags;
    reset_class event;
};

constexpr std::array<named_flag, 10> kFlags{{{RCC_RSR_CPURSTF, "cpu"},
                                             {RCC_RSR_D1RSTF, "d1"},
                                             {RCC_RSR_D2RSTF, "d2"},
                                             {RCC_RSR_BORRSTF, "bor"},
                                             {RCC_RSR_PINRSTF, "pin"},
                                             {RCC_RSR_PORRSTF, "por"},
                                             {RCC_RSR_SFTRSTF, "sft"},
                                             {RCC_RSR_IWDG1RSTF, "iwdg1"},
                                             {RCC_RSR_WWDG1RSTF, "wwdg1"},
                                             {RCC_RSR_LPWRRSTF, "lpwr"}}};

constexpr std::uint32_t kPinCpu = RCC_RSR_PINRSTF | RCC_RSR_CPURSTF;

// RM0433 section 8.4.4, "Reset source identification": the flags each reset
// source leaves set in RCC_RSR.
constexpr std::array<event_signature, 11> kEvents{{{RCC_RSR_PORRSTF | RCC_RSR_BORRSTF | RCC_RSR_D2RSTF | RCC_RSR_D1RSTF | kPinCpu, reset_class::power_on},
                                                   {kPinCpu, reset_class::pin},
                                                   {RCC_RSR_BORRSTF | kPinCpu, reset_class::brownout},
                                                   {RCC_RSR_SFTRSTF | kPinCpu, reset_class::software},
                                                   {RCC_RSR_CPURSTF, reset_class::cpu},
                                                   {RCC_RSR_WWDG1RSTF, reset_class::window_watchdog},
                                                   {RCC_RSR_WWDG1RSTF | kPinCpu, reset_class::window_watchdog},
                                                   {RCC_RSR_IWDG1RSTF | kPinCpu, reset_class::independent_watchdog},
                                                   {RCC_RSR_D1RSTF, reset_class::d1_standby_exit},
                                                   {RCC_RSR_D2RSTF, reset_class::d2_standby_exit},
                                                   {RCC_RSR_LPWRRSTF | kPinCpu, reset_class::low_power}}};

constexpr std::uint32_t known_flags() noexcept
{
    std::uint32_t mask = 0;
    for(const named_flag &flag : kFlags)
        mask |= flag.mask;
    return mask;
}

reset_class classify(std::uint32_t flags) noexcept
{
    for(const event_signature &signature : kEvents)
        if(signature.flags == flags)
            return signature.event;
    return reset_class::unmatched;
}

const char *name_of(reset_class event) noexcept
{
    constexpr std::array<const char *, 11> names{"power_on",        "pin",       "brownout", "software", "cpu", "window_watchdog", "independent_watchdog", "d1_standby_exit",
                                                 "d2_standby_exit", "low_power", "unmatched"};
    return names[static_cast<std::size_t>(event)];
}

void print_flags(std::uint32_t status)
{
    const char *separator = "";
    for(const named_flag &flag : kFlags)
    {
        if((status & flag.mask) == 0U)
            continue;
        std::printf("%s%s", separator, flag.name);
        separator = "+";
    }
    if(*separator == '\0')
        std::printf("none");
}

}

reset_cause read_reset_cause() noexcept
{
    const std::uint32_t status = RCC->RSR;
    return {status, classify(status & known_flags())};
}

// RM0433 describes RMVF as set and reset by software, so it is written back to
// zero rather than left set holding the flags cleared through the next reset.
void clear_reset_cause() noexcept
{
    RCC->RSR |= RCC_RSR_RMVF;
    RCC->RSR &= ~RCC_RSR_RMVF;
}

void print_reset_cause(const reset_cause &cause)
{
    std::printf("reset=%s reset_flags=", name_of(cause.event));
    print_flags(cause.status);
    std::printf(" reset_rsr=0x%08" PRIX32, cause.status);
}

}
