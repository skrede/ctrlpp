#include "posture.h"

#include "stm32h7xx.h"

#include <cstdio>
#include <cstdint>
#include <cinttypes>

// The heap's base, placed by the linker script.
extern "C" char end[];

namespace {

constexpr std::uint64_t kCounterStates = std::uint64_t{1} << 32;
constexpr std::uint64_t kMillisecond   = 1000;

const char *clock_source() noexcept
{
    switch(RCC->CFGR & RCC_CFGR_SWS)
    {
        case RCC_CFGR_SWS_HSI:
            return "hsi";
        case RCC_CFGR_SWS_CSI:
            return "csi";
        case RCC_CFGR_SWS_HSE:
            return "hse";
        case RCC_CFGR_SWS_PLL1:
            return "pll1";
    }
    return "unknown";
}

// The device header gives each region's base but not its size; a valid address
// lies below the next region's base.
const char *region_of(std::uintptr_t address) noexcept
{
    if(address >= D1_AXIFLASH_BASE && address <= FLASH_END)
        return "flash";
    if(address >= D1_DTCMRAM_BASE && address < D1_AXISRAM_BASE)
        return "dtcm";
    if(address >= D1_AXISRAM_BASE && address < D2_AHBSRAM_BASE)
        return "axi_sram";
    return "other";
}

template<class T>
const char *region_of(T *object) noexcept
{
    return region_of(reinterpret_cast<std::uintptr_t>(object));
}

const char *bit_state(std::uint32_t reg, std::uint32_t mask) noexcept
{
    return (reg & mask) != 0U ? "on" : "off";
}

const char *optimization() noexcept
{
#if defined(__OPTIMIZE_SIZE__)
    return "size";
#elif defined(__OPTIMIZE__)
    return "speed";
#else
    return "none";
#endif
}

}

namespace ctrlpp {

void report_posture()
{
    SystemCoreClockUpdate();
    const std::uint32_t ccr      = SCB->CCR;
    const std::uint32_t latency  = (FLASH->ACR & FLASH_ACR_LATENCY_Msk) >> FLASH_ACR_LATENCY_Pos;
    const std::uint64_t clock_hz = SystemCoreClock != 0U ? SystemCoreClock : 1U;
    const std::uint32_t wrap_ms  = static_cast<std::uint32_t>(kCounterStates * kMillisecond / clock_hz);
    char frame_marker            = 0;
    std::printf("[posture] core_hz=%" PRIu32 " axi_hz=%" PRIu32 " clock_source=%s icache=%s dcache=%s flash_latency_ws=%" PRIu32 " cycle_wrap_ms=%" PRIu32, SystemCoreClock,
                SystemD2Clock, clock_source(), bit_state(ccr, SCB_CCR_IC_Msk), bit_state(ccr, SCB_CCR_DC_Msk), latency, wrap_ms);
    std::printf(" code=%s data=%s stack=%s heap=%s optimize=%s compiler=gcc-%s\n", region_of(&report_posture), region_of(&SystemCoreClock), region_of(&frame_marker), region_of(end),
                optimization(), __VERSION__);
}

}
