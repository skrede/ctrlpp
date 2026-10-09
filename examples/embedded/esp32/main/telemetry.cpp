#include "telemetry.h"

#include "golden_reference.h"

#include "sdkconfig.h"
#include "esp_rom_sys.h"
#include "esp_app_desc.h"
#include "esp_idf_version.h"

#include "driver/uart.h"

#include <span>
#include <array>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <cinttypes>
#include <string_view>

namespace ctrlpp {

namespace {

// Telemetry UART (UART2): ESP TX GPIO5 -> USB RX, ESP RX GPIO4 <- USB TX. The
// Olimex ESP32-POE boards drive the Ethernet PHY clock from GPIO17, and GPIO3 is
// the console UART's RX from the on-board USB bridge, so neither may carry it.
// GPIO5 is a strapping pin that must be high at boot, which an idle UART line is.
constexpr uart_port_t kTelemetryUart = UART_NUM_2;
constexpr std::int32_t kTxPin        = 5;
constexpr std::int32_t kRxPin        = 4;
constexpr std::int32_t kBaud         = 115200;

constexpr std::size_t kShaBytes = sizeof(esp_app_desc_t::app_elf_sha256);

// The watchpoint makes up to 60 bytes of every task's stack unusable, so a stack
// figure is only comparable with one taken under the same setting.
#ifdef CONFIG_FREERTOS_WATCHPOINT_END_OF_STACK
constexpr bool kStackWatchpoint = true;
#else
constexpr bool kStackWatchpoint = false;
#endif

void write(std::string_view text)
{
    uart_write_bytes(kTelemetryUart, text.data(), text.size());
}

std::size_t append(std::span<char> line, std::size_t used, std::string_view text)
{
    const std::size_t n = std::min(text.size(), line.size() - 1 - used);
    std::copy_n(text.data(), n, line.data() + used);
    line[used + n] = '\0';
    return used + n;
}

// The sha256 digest of the ELF the image was made from, which the build writes into
// the image's descriptor; hashing the built ELF on the host reproduces it.
std::array<char, 2 * kShaBytes + 1> elf_sha256_text()
{
    constexpr std::string_view digits = "0123456789abcdef";
    const esp_app_desc_t *desc        = esp_app_get_description();
    std::array<char, 2 * kShaBytes + 1> text{};
    for(std::size_t i = 0; i < kShaBytes; ++i)
    {
        text[2 * i]     = digits[desc->app_elf_sha256[i] >> 4U];
        text[2 * i + 1] = digits[desc->app_elf_sha256[i] & 0xFU];
    }
    return text;
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

void telemetry_init()
{
    uart_config_t cfg{};
    cfg.baud_rate  = kBaud;
    cfg.data_bits  = UART_DATA_8_BITS;
    cfg.parity     = UART_PARITY_DISABLE;
    cfg.stop_bits  = UART_STOP_BITS_1;
    cfg.flow_ctrl  = UART_HW_FLOWCTRL_DISABLE;
    cfg.source_clk = UART_SCLK_DEFAULT;

    uart_driver_install(kTelemetryUart, 512, 0, 0, nullptr, 0);
    uart_param_config(kTelemetryUart, &cfg);
    uart_set_pin(kTelemetryUart, kTxPin, kRxPin, UART_PIN_NO_CHANGE, UART_PIN_NO_CHANGE);
}

void send_line(stream_tag tag, const char *body)
{
    write(tag_text(tag));
    write(body);
    write("\n");
}

void emit_identity()
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "build_id=%s", elf_sha256_text().data());
    send_line(stream_tag::meta, body);
    send_line(stream_tag::meta, "kind=uninstrumented");
    std::snprintf(body, sizeof(body), "board=ESP32-POE-EA target=%s loop=freertos-task-50hz scalar=float", CONFIG_IDF_TARGET);
    send_line(stream_tag::meta, body);
}

void emit_posture(const reset_cause &cause)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "posture reset=%s reset_reason=%" PRId32 " reset_rom=%" PRIu32 " core_mhz=%" PRIu32 " stack_watchpoint=%s optimize=%s compiler=gcc-%s idf=%s",
                  name_of(cause.event), cause.reason, cause.rom_code, esp_rom_get_cpu_ticks_per_us(), kStackWatchpoint ? "on" : "off", optimization(), __VERSION__,
                  esp_get_idf_version());
    send_line(stream_tag::meta, body);
}

void emit_canaries()
{
    send_line(stream_tag::canary, "allocation UNAVAILABLE -- this image carries no allocation sensor; no allocation figure is printed");
    send_line(stream_tag::canary, "cycle_counter UNAVAILABLE -- this image carries no cycle counter; no timing figure is printed");
}

void emit_header(claim_class claim)
{
    char body[kLineBytes] = {};
    std::size_t used      = 0;
    std::string_view separator;
    for(const claim_field &field : kClaimFields)
    {
        if(field.claim != claim)
            continue;
        used      = append(body, used, separator);
        used      = append(body, used, field.name);
        separator = ",";
    }
    send_line(tag_of(claim), body);
}

void emit_exact_row(std::int32_t k, const std::array<float, kExactValues> &values)
{
    static_assert(kExactValues == 4, "the row's format names four values");
    std::array<std::array<char, kBitsTextLength<float>>, kExactValues> bits{};
    for(std::size_t i = 0; i < kExactValues; ++i)
        render_bits(values[i], bits[i]);
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "%" PRId32 ",%s,%s,%s,%s", k, bits[0].data(), bits[1].data(), bits[2].data(), bits[3].data());
    send_line(stream_tag::exact, body);
}

void emit_gain(float k0, float k1)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "gain device=%.17g,%.17g host=%.17g,%.17g", static_cast<double>(k0), static_cast<double>(k1), kHostK0, kHostK1);
    send_line(stream_tag::verdict, body);
}

void emit_stack(const char *region, std::uint32_t free_bytes, std::uint32_t stack_bytes)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "stack region=%s used=%" PRIu32 " free=%" PRIu32 " of=%" PRIu32 " bytes", region, stack_bytes - free_bytes, free_bytes, stack_bytes);
    send_line(stream_tag::verdict, body);
}

}
