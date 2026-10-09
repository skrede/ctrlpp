#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_CLOCK_TREE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_CLOCK_TREE_H

#include <cstdint>

namespace ctrlpp {

// The USART234578SEL code for hsi_ker_ck (RM0433, RCC_D2CCIP2R).
inline constexpr std::uint32_t kUsartKernelHsi = 3U;

// DBGMCU_IDC REV_ID (RM0433, DBGMCU chapter): 0x1003 is revision Y, 0x2003
// revision V.
std::uint32_t silicon_revision() noexcept;

// Appends the fields the tuned posture adds to the posture line: the silicon
// revision that chose the clock, the voltage scale, the PLL's reference and the
// console's kernel clock, each read back from the hardware.
void print_clock_tree();

}

#endif
