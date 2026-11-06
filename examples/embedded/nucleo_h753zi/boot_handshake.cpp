#include "telemetry.h"
#include "boot_handshake.h"

#include "stm32h7xx.h"

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace
{

struct handshake_outcome
{
    bool received;
    std::uint32_t waited_ms;
};

// PD9 is USART3_RX, bridged to the ST-LINK's virtual COM port beside the TX pin.
constexpr std::uint32_t rx_pin = 9;
constexpr std::uint32_t usart3_function = 7;
constexpr std::uint32_t milliseconds_per_s = 1000;
constexpr std::uint32_t damaged = USART_ISR_FE | USART_ISR_NE | USART_ISR_PE;
constexpr std::uint32_t clear_errors = USART_ICR_ORECF | USART_ICR_FECF | USART_ICR_NECF | USART_ICR_PECF;

// Static, so main's frame, which lies beneath every stack figure, does not grow.
handshake_outcome outcome = {false, 0};

void receiver_on() noexcept
{
    GPIOD->MODER = (GPIOD->MODER & ~(0x3u << (rx_pin * 2))) | (0x2u << (rx_pin * 2));
    GPIOD->AFR[1] = (GPIOD->AFR[1] & ~(0xFu << ((rx_pin - 8) * 4))) | (usart3_function << ((rx_pin - 8) * 4));
    USART3->ICR = clear_errors;
    USART3->CR1 |= USART_CR1_RE;
}

// A byte damaged in transit does not count, so line noise while the ST-LINK
// powers up cannot stand in for the host. An overrun stalls the receiver until
// it is cleared, so every pass clears it.
bool handshake_arrived() noexcept
{
    const std::uint32_t status = USART3->ISR;
    USART3->ICR = clear_errors;
    if((status & USART_ISR_RXNE_RXFNE) == 0U)
        return false;
    const auto byte = static_cast<std::uint8_t>(USART3->RDR);
    return (status & damaged) == 0U && byte == ctrlpp::boot_handshake_byte;
}

}

namespace ctrlpp
{

// SysTick counts milliseconds with its interrupt left disabled, polled through
// its count flag, so the wait needs no interrupt and does not touch the cycle
// counter the timing image proves afterwards.
void await_boot_handshake() noexcept
{
    SystemCoreClockUpdate();
    receiver_on();
    SysTick->LOAD = SystemCoreClock / milliseconds_per_s - 1U;
    SysTick->VAL = 0U;
    SysTick->CTRL = SysTick_CTRL_CLKSOURCE_Msk | SysTick_CTRL_ENABLE_Msk;
    bool received = handshake_arrived();
    std::uint32_t ms = 0;
    while(!received && ms < boot_handshake_timeout_ms)
    {
        if((SysTick->CTRL & SysTick_CTRL_COUNTFLAG_Msk) != 0U)
            ++ms;
        received = handshake_arrived();
    }
    SysTick->CTRL = 0U;
    USART3->CR1 &= ~USART_CR1_RE;
    outcome = {received, ms};
}

void emit_boot_handshake()
{
    open_line(stream_tag::meta);
    std::printf("boot_handshake state=%s waited_ms=%" PRIu32 " timeout_ms=%" PRIu32 "\n", outcome.received ? "received" : "timeout", outcome.waited_ms, boot_handshake_timeout_ms);
}

}
