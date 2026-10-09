#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_RESET_REASON_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_RESET_REASON_H

#include <cstdint>

namespace ctrlpp {

enum class reset_class : std::uint8_t
{
    power_on,
    pin,
    software,
    brownout,
    panic,
    interrupt_watchdog,
    task_watchdog,
    watchdog,
    deep_sleep_exit,
    sdio,
    usb,
    jtag,
    efuse,
    power_glitch,
    cpu_lockup,
    unknown,
    unmatched
};

struct reset_cause
{
    std::int32_t reason;
    std::uint32_t rom_code;
    reset_class event;
};

// Nothing is cleared after the read, unlike a latched flag register: the
// hardware code is rewritten by every reset, and the framework clears its own
// hint while it reads it at startup.
reset_cause read_reset_cause() noexcept;

const char *name_of(reset_class event) noexcept;

}

#endif
