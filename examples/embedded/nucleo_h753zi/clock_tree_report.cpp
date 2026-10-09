#include "clock_tree.h"

#include "stm32h7xx.h"

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace {

std::uint32_t voltage_scale() noexcept
{
    const bool syscfg_clocked = (RCC->APB4ENR & RCC_APB4ENR_SYSCFGEN) != 0U;
    if(syscfg_clocked && (SYSCFG->PWRCR & SYSCFG_PWRCR_ODEN) != 0U)
        return 0U;
    return 4U - ((PWR->D3CR & PWR_D3CR_VOS) >> PWR_D3CR_VOS_Pos);
}

const char *pll_source() noexcept
{
    switch(RCC->PLLCKSELR & RCC_PLLCKSELR_PLLSRC)
    {
        case RCC_PLLCKSELR_PLLSRC_HSI:
            return "hsi";
        case RCC_PLLCKSELR_PLLSRC_CSI:
            return "csi";
        case RCC_PLLCKSELR_PLLSRC_HSE:
            return "hse";
    }
    return "none";
}

const char *console_kernel() noexcept
{
    const std::uint32_t selected = (RCC->D2CCIP2R & RCC_D2CCIP2R_USART28SEL) >> RCC_D2CCIP2R_USART28SEL_Pos;
    return selected == 0U ? "pclk1" : selected == ctrlpp::kUsartKernelHsi ? "hsi" : "other";
}

}

namespace ctrlpp {

void print_clock_tree()
{
    std::printf(" silicon_rev=0x%04" PRIX32 " vos=%" PRIu32 " pll_source=%s usart3_kernel=%s", silicon_revision(), voltage_scale(), pll_source(), console_kernel());
}

}
