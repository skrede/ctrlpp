#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_BOOT_HANDSHAKE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_BOOT_HANDSHAKE_H

#include <cstdint>

namespace ctrlpp {

// ASCII SYN: a control byte that the AT probes a modem manager writes to a
// newly attached serial port never contain, so a probe cannot release the wait.
inline constexpr std::uint8_t kBootHandshakeByte = 0x16;

// The wait must outlast the host's side of a restart: after drag-and-drop
// programming the ST-LINK took about 11 s to re-enumerate before the capture
// could reopen its port. A longer wait costs nothing with the host attached,
// whose byte lands about 1 s after the port returns, so the default leaves a
// wide margin. Configured with CTRLPP_MCU_BOOT_HANDSHAKE_MS.
inline constexpr std::uint32_t kBootHandshakeTimeoutMs = CTRLPP_MCU_BOOT_HANDSHAKE_MS;

// Waits for the host's byte on USART3, then turns the receiver off again. The
// console must already be up. Runs before the cycle counter is enabled and
// leaves the core's interrupts off, so no figure in the report contains it.
void await_boot_handshake() noexcept;

// A run whose wait expired began reporting with no host attached, so its head
// may be lost; the line says which happened and how long the wait took.
void emit_boot_handshake();

}

#endif
