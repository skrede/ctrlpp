#include "reset_reason.h"

#include "esp_system.h"
#include "esp_rom_sys.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

namespace {

struct reason_class
{
    esp_reset_reason_t reason;
    reset_class event;
};

// The class is the framework's, not the reset register's: the register reports a
// panic, an esp_restart and a task-watchdog reset alike as a software reset of
// the core, and only the hint the framework keeps across the reset tells them
// apart. On this part the chip-enable pin and a power cycle latch the same
// power-on code, so power_on covers both, and the framework documents
// ESP_RST_EXT as not applicable here.
constexpr std::array<reason_class, 16> kClasses{{{ESP_RST_POWERON, reset_class::power_on},
                                                 {ESP_RST_EXT, reset_class::pin},
                                                 {ESP_RST_SW, reset_class::software},
                                                 {ESP_RST_BROWNOUT, reset_class::brownout},
                                                 {ESP_RST_PANIC, reset_class::panic},
                                                 {ESP_RST_INT_WDT, reset_class::interrupt_watchdog},
                                                 {ESP_RST_TASK_WDT, reset_class::task_watchdog},
                                                 {ESP_RST_WDT, reset_class::watchdog},
                                                 {ESP_RST_DEEPSLEEP, reset_class::deep_sleep_exit},
                                                 {ESP_RST_SDIO, reset_class::sdio},
                                                 {ESP_RST_USB, reset_class::usb},
                                                 {ESP_RST_JTAG, reset_class::jtag},
                                                 {ESP_RST_EFUSE, reset_class::efuse},
                                                 {ESP_RST_PWR_GLITCH, reset_class::power_glitch},
                                                 {ESP_RST_CPU_LOCKUP, reset_class::cpu_lockup},
                                                 {ESP_RST_UNKNOWN, reset_class::unknown}}};

reset_class classify(esp_reset_reason_t reason) noexcept
{
    for(const reason_class &entry : kClasses)
        if(entry.reason == reason)
            return entry.event;
    return reset_class::unmatched;
}

}

reset_cause read_reset_cause() noexcept
{
    const esp_reset_reason_t reason = esp_reset_reason();
    return {static_cast<std::int32_t>(reason), static_cast<std::uint32_t>(esp_rom_get_reset_reason(0)), classify(reason)};
}

const char *name_of(reset_class event) noexcept
{
    constexpr std::array<const char *, 17> names{"power_on", "pin", "software", "brownout", "panic",        "interrupt_watchdog", "task_watchdog", "watchdog", "deep_sleep_exit",
                                                 "sdio",     "usb", "jtag",     "efuse",    "power_glitch", "cpu_lockup",         "unknown",       "unmatched"};
    return names[static_cast<std::size_t>(event)];
}

}
