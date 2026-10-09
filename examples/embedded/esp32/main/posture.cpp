#include "posture.h"
#include "telemetry.h"
#include "cycle_counter.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_cpu.h"
#include "sdkconfig.h"
#include "esp_rom_sys.h"
#include "esp_task_wdt.h"
#include "esp_idf_version.h"
#include "esp_memory_utils.h"

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace {

// The watchpoint makes up to 60 bytes of every task's stack unusable, so a stack
// figure is only comparable with one taken under the same setting.
#ifdef CONFIG_FREERTOS_WATCHPOINT_END_OF_STACK
constexpr bool kStackWatchpoint = true;
#else
constexpr bool kStackWatchpoint = false;
#endif

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

// The framework's default build keeps C assertions, and with them Eigen's own
// checks; the timing image is built without them, the evidence image with them,
// and a figure from one build is not a figure of the other.
#ifdef NDEBUG
constexpr bool kAssertions = false;
#else
constexpr bool kAssertions = true;
#endif

const char *scheduler_state() noexcept
{
    switch(xTaskGetSchedulerState())
    {
        case taskSCHEDULER_RUNNING:
            return "running";
        case taskSCHEDULER_SUSPENDED:
            return "suspended";
        case taskSCHEDULER_NOT_STARTED:
            return "not_started";
    }
    return "unknown";
}

// Read from the framework rather than the configuration, since the timing image
// stops the watchdog the configuration starts.
const char *task_watchdog() noexcept
{
    return esp_task_wdt_status(nullptr) == ESP_ERR_INVALID_STATE ? "off" : "on";
}

// Code outside instruction RAM is fetched from flash through the cache, so its
// first execution pays for the misses and later ones may not.
const char *code_region() noexcept
{
    const void *code = reinterpret_cast<const void *>(&ctrlpp::emit_posture);
    if(esp_ptr_in_iram(code))
        return "iram";
    return esp_ptr_executable(code) ? "flash" : "other";
}

}

namespace ctrlpp {

void emit_posture(const reset_cause &cause)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body),
                  "posture reset=%s reset_reason=%" PRId32 " reset_rom=%" PRIu32 " core_mhz=%" PRIu32 " cycle_bits=%" PRIu32 " cycle_wrap_ms=%" PRIu32
                  " stack_watchpoint=%s assertions=%s optimize=%s compiler=gcc-%s idf=%s",
                  name_of(cause.event), cause.reason, cause.rom_code, esp_rom_get_cpu_ticks_per_us(), kCounterBits, cycle_wrap_ms(), kStackWatchpoint ? "on" : "off",
                  kAssertions ? "on" : "off", optimization(), __VERSION__, esp_get_idf_version());
    send_line(stream_tag::meta, body);
    std::snprintf(body, sizeof(body),
                  "schedule tick_hz=%" PRIu32 " task_core=%" PRId32 " task_priority=%" PRIu32 " scheduler=%s intr_enabled=0x%08" PRIx32
                  " task_wdt=%s code=%s flash_mode=%s flash_freq=%s",
                  static_cast<std::uint32_t>(configTICK_RATE_HZ), static_cast<std::int32_t>(xPortGetCoreID()), static_cast<std::uint32_t>(uxTaskPriorityGet(nullptr)), scheduler_state(),
                  esp_cpu_intr_get_enabled_mask(), task_watchdog(), code_region(), CONFIG_ESPTOOLPY_FLASHMODE, CONFIG_ESPTOOLPY_FLASHFREQ);
    send_line(stream_tag::meta, body);
}

}
