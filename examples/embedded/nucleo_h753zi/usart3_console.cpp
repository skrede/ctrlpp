// Minimal USART3 console for headless CSV capture over the ST-Link VCP.
//
// The NUCLEO-144 ST-Link virtual COM port is wired to USART3 (PD8=TX, PD9=RX)
// by the on-board solder bridges, so routing the console descriptors' _write
// (syscalls.cpp) out USART3 makes std::printf reach the host's /dev/ttyACM0 with
// no debugger/semihosting. After the CMSIS SystemInit (which does not configure
// the PLL) the H753 runs on HSI = 64 MHz, so PCLK1 = 64 MHz and USART3's kernel
// clock = 64 MHz; BRR = 64e6 / 115200 ~= 556 for oversampling-by-16. The tuned
// posture keeps that kernel clock by moving USART3 onto HSI (clock_tree.cpp).

#include "stm32h7xx.h"

#include <span>
#include <cstddef>

namespace ctrlpp
{

void usart3_console_init() noexcept
{
    // Read the enable register back after setting each clock-enable bit so the
    // enable has taken effect before the peripheral registers are touched.
    RCC->AHB4ENR |= RCC_AHB4ENR_GPIODEN;
    (void)RCC->AHB4ENR;
    RCC->APB1LENR |= RCC_APB1LENR_USART3EN;
    (void)RCC->APB1LENR;

    // PD8 = USART3_TX, alternate function 7. PD9, the receiver, is configured
    // only for the boot handshake, by boot_handshake.cpp.
    GPIOD->MODER = (GPIOD->MODER & ~(0x3u << (8 * 2))) | (0x2u << (8 * 2));
    GPIOD->AFR[1] = (GPIOD->AFR[1] & ~(0xFu << ((8 - 8) * 4)))
                  | (0x7u << ((8 - 8) * 4));

    USART3->CR1 = 0;
    USART3->BRR = 556;
    USART3->CR1 = USART_CR1_TE | USART_CR1_UE;
}

inline void usart3_putc(char c) noexcept
{
    while(!(USART3->ISR & USART_ISR_TXE_TXFNF)) { }
    USART3->TDR = static_cast<uint8_t>(c);
}

void usart3_write(std::span<const char> text) noexcept
{
    for(const char c : text)
    {
        if(c == '\n')
            usart3_putc('\r');
        usart3_putc(c);
    }
}

}
