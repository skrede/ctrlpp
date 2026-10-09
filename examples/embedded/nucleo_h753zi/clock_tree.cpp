#include "clock_tree.h"

#include "stm32h7xx.h"

#include <cstdint>
#include <optional>

extern "C" void __real_SystemInit();

namespace {

struct clock_plan
{
    std::uint32_t pll_n;
    std::uint32_t flash_latency;
    std::uint32_t flash_write_delay;
    bool overdrive;
};

// Revision Y reaches 400 MHz in VOS1, revision V 480 MHz in VOS0, the overdrive
// scale (DS12117, general operating conditions). Any other revision keeps the
// reset clock tree.
constexpr std::uint32_t kRevisionY = 0x1003U;
constexpr std::uint32_t kRevisionV = 0x2003U;

// The core is VCO / 2 with VCO = 16 MHz * N, inside the wide range 192 to
// 960 MHz; the AXI clock is the core's half, 240 or 200 MHz. Flash rows from
// RM0433 Table 17, "FLASH recommended number of wait states and programming
// delay": 4 WS and WRHIGHFREQ 2 up to 240 MHz in VOS0, 2 WS and 2 up to 210 MHz
// in VOS1.
constexpr clock_plan kRevisionVPlan = {60U, 4U, 2U, true};
constexpr clock_plan kRevisionYPlan = {50U, 2U, 2U, false};

constexpr std::uint32_t kPllM          = 4U;
constexpr std::uint32_t kPllP          = 2U;
constexpr std::uint32_t kInput8To16    = 3U;
constexpr std::uint32_t kVoltageScale1 = 3U;
constexpr std::uint32_t kSupplyConfig  = PWR_CR3_SCUEN | PWR_CR3_LDOEN | PWR_CR3_BYPASS;
constexpr std::uint32_t kPllShape      = RCC_PLLCFGR_PLL1FRACEN | RCC_PLLCFGR_PLL1VCOSEL | RCC_PLLCFGR_PLL1RGE | RCC_PLLCFGR_DIVQ1EN | RCC_PLLCFGR_DIVR1EN;
constexpr std::uint32_t kFlashTiming   = FLASH_ACR_LATENCY | FLASH_ACR_WRHIGHFREQ;

std::optional<clock_plan> plan_for(std::uint32_t revision) noexcept
{
    if(revision == kRevisionV)
        return kRevisionVPlan;
    if(revision == kRevisionY)
        return kRevisionYPlan;
    return std::nullopt;
}

// The supply configuration is writable once per power-on (RM0433, PWR_CR3
// SCUEN); after any other reset it holds what the previous run locked.
bool lock_ldo_supply() noexcept
{
    if((PWR->CR3 & PWR_CR3_SCUEN) == 0U)
        return (PWR->CR3 & kSupplyConfig) == PWR_CR3_LDOEN;
    PWR->CR3 = (PWR->CR3 & ~kSupplyConfig) | PWR_CR3_LDOEN;
    while((PWR->CSR1 & PWR_CSR1_ACTVOSRDY) == 0U)
    {
    }
    return true;
}

void await_voltage() noexcept
{
    while((PWR->D3CR & PWR_D3CR_VOSRDY) == 0U)
    {
    }
}

// VOS0 is entered from VOS1 by the SYSCFG overdrive bit, the sequence of
// __HAL_PWR_VOLTAGESCALING_CONFIG in STM32CubeH7's stm32h7xx_hal_pwr.h.
void raise_voltage(bool overdrive) noexcept
{
    PWR->D3CR = (PWR->D3CR & ~PWR_D3CR_VOS) | (kVoltageScale1 << PWR_D3CR_VOS_Pos);
    await_voltage();
    if(!overdrive)
        return;
    RCC->APB4ENR |= RCC_APB4ENR_SYSCFGEN;
    (void)RCC->APB4ENR;
    SYSCFG->PWRCR |= SYSCFG_PWRCR_ODEN;
    (void)SYSCFG->PWRCR;
    await_voltage();
}

void start_pll1(std::uint32_t n) noexcept
{
    RCC->PLLCKSELR = (RCC->PLLCKSELR & ~(RCC_PLLCKSELR_PLLSRC | RCC_PLLCKSELR_DIVM1)) | RCC_PLLCKSELR_PLLSRC_HSI | (kPllM << RCC_PLLCKSELR_DIVM1_Pos);
    RCC->PLLCFGR   = (RCC->PLLCFGR & ~kPllShape) | (kInput8To16 << RCC_PLLCFGR_PLL1RGE_Pos) | RCC_PLLCFGR_DIVP1EN;
    RCC->PLL1DIVR  = (RCC->PLL1DIVR & ~(RCC_PLL1DIVR_N1 | RCC_PLL1DIVR_P1)) | ((n - 1U) << RCC_PLL1DIVR_N1_Pos) | ((kPllP - 1U) << RCC_PLL1DIVR_P1_Pos);
    RCC->CR |= RCC_CR_PLL1ON;
    while((RCC->CR & RCC_CR_PLL1RDY) == 0U)
    {
    }
}

// The AXI and AHB buses run at half the core and every APB at a quarter, which
// is each bus's ceiling under either plan: 240 and 120 MHz in VOS0, 200 and
// 100 MHz in VOS1.
void set_bus_prescalers() noexcept
{
    RCC->D1CFGR = RCC_D1CFGR_HPRE_DIV2 | RCC_D1CFGR_D1PPRE_DIV2;
    RCC->D2CFGR = RCC_D2CFGR_D2PPRE1_DIV2 | RCC_D2CFGR_D2PPRE2_DIV2;
    RCC->D3CFGR = RCC_D3CFGR_D3PPRE_DIV2;
}

void set_flash_wait_states(const clock_plan &plan) noexcept
{
    FLASH->ACR = (FLASH->ACR & ~kFlashTiming) | plan.flash_latency | (plan.flash_write_delay << FLASH_ACR_WRHIGHFREQ_Pos);
    while((FLASH->ACR & FLASH_ACR_LATENCY) != plan.flash_latency)
    {
    }
}

void switch_to_pll1() noexcept
{
    RCC->CFGR = (RCC->CFGR & ~RCC_CFGR_SW) | RCC_CFGR_SW_PLL1;
    while((RCC->CFGR & RCC_CFGR_SWS) != RCC_CFGR_SWS_PLL1)
    {
    }
}

// usart3_console.cpp divides a 64 MHz kernel clock, which APB1 no longer is.
void keep_console_on_hsi() noexcept
{
    RCC->D2CCIP2R = (RCC->D2CCIP2R & ~RCC_D2CCIP2R_USART28SEL) | (ctrlpp::kUsartKernelHsi << RCC_D2CCIP2R_USART28SEL_Pos);
}

void raise_core_clock(const clock_plan &plan) noexcept
{
    raise_voltage(plan.overdrive);
    start_pll1(plan.pll_n);
    set_bus_prescalers();
    set_flash_wait_states(plan);
    switch_to_pll1();
    keep_console_on_hsi();
}

// TCM accesses bypass the L1 caches (Arm Cortex-M7 TRM, memory system), so the
// stack, .data, .bss and the sensor's counters in DTCM are unaffected. Only the
// core touches the cached AXI SRAM heap, so it needs no maintenance.
void enable_caches() noexcept
{
    SCB_EnableICache();
    SCB_EnableDCache();
}

}

// Runs from the reset handler before .data is copied and .bss zeroed, so it
// reads and writes nothing but registers and its own stack.
extern "C" void __wrap_SystemInit()
{
    __real_SystemInit();
    const std::optional<clock_plan> plan = plan_for(ctrlpp::silicon_revision());
    if(plan.has_value() && lock_ldo_supply())
        raise_core_clock(*plan);
    enable_caches();
}

namespace ctrlpp {

std::uint32_t silicon_revision() noexcept
{
    return (DBGMCU->IDCODE & DBGMCU_IDCODE_REV_ID_Msk) >> DBGMCU_IDCODE_REV_ID_Pos;
}

}
