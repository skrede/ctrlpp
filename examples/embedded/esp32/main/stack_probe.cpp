// Measures on the target how large a state dimension the discrete solve fits in a task's stack, where a
// host compiler's -fstack-usage report only predicts it for another architecture and calling convention.
// Pass one gives every dimension a generous stack, so each reports what it needs; pass two gives each the
// control task's stack size, so the first dimension that does not fit establishes the supported maximum.
// Each dimension runs in its own task, since the high-water mark is a running minimum over a task's life,
// and in increasing order, each reporting before the next starts, so an overflow leaves a complete record
// of every dimension below it.

#include "stack_probe.h"
#include "stack_probe_dimension.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "freertos/semphr.h"

#include "esp_log.h"
#include "esp_rom_sys.h"

#include <cstdint>
#include <cinttypes>

namespace
{

using ctrlpp::probe::probe_tag;

// Deliberately generous, so that pass one measures what each dimension NEEDS
// rather than only whether it fitted the control task's stack. A dimension that
// overflows here would be reported the same way, by name.
constexpr std::uint32_t headroom_stack_bytes = 48 * 1024;

void sweep(std::uint32_t stack, SemaphoreHandle_t done)
{
    using ctrlpp::probe::run_one;
    run_one<2, 1>(stack, done);
    run_one<3, 1>(stack, done);
    run_one<4, 2>(stack, done);
    run_one<5, 1>(stack, done);
    run_one<6, 2>(stack, done);
    run_one<8, 2>(stack, done);
}

void headroom_pass(SemaphoreHandle_t done)
{
    ESP_LOGI(probe_tag, "PASS 1 -- requirement per dimension, %" PRIu32 " byte stack, scalar = float", headroom_stack_bytes);
    sweep(headroom_stack_bytes, done);
    ESP_LOGI(probe_tag, "PASS 1 complete -- every dimension returned");
}

// Expected to end in a named overflow, and that is the point: it establishes
// the supported maximum by reaching the first dimension that does not fit.
void limit_pass(SemaphoreHandle_t done)
{
    ESP_LOGI(probe_tag, "PASS 2 -- supported maximum at the control task's own %" PRIu32 " byte stack", ctrlpp::control_task_stack_bytes);
    sweep(ctrlpp::control_task_stack_bytes, done);
    ESP_LOGI(probe_tag, "PASS 2 complete -- no dimension overflowed, which was not expected");
}

void stack_probe_task(void*)
{
    SemaphoreHandle_t done = xSemaphoreCreateBinary();
    headroom_pass(done);
    limit_pass(done);
    vTaskDelete(nullptr);
}

}

namespace ctrlpp
{

// The driver task itself needs almost no stack: every solve happens in a task
// this one creates, and it only waits on them.
void start_stack_probe()
{
    xTaskCreatePinnedToCore(stack_probe_task, "ctrlpp_probe", 3 * 1024, nullptr, 4, nullptr, 0);
}

}

// Overrides the framework's weak default so an overflow names its task rather than resetting the board
// unexplained. esp_rom_printf writes the console UART directly, since the log path takes locks and
// allocates and a task that has just overrun its stack cannot be trusted to complete it; its conversions
// are the ROM printer's own, hence %u.
extern "C" void vApplicationStackOverflowHook(TaskHandle_t, char* name)
{
    esp_rom_printf("STACK OVERFLOW in task '%s' -- the dimension above the last "
                   "reported probe line does not fit %u bytes\n",
            name != nullptr ? name : "<unnamed>", static_cast<unsigned>(ctrlpp::control_task_stack_bytes));
}
