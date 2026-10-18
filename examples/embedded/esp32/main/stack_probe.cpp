// How large a state dimension the discrete solve fits in a task's stack.
// Every stack figure this project publishes for the Riccati path so far comes
// from a host compiler's -fstack-usage report. That is a PREDICTION about a
// different architecture, a different calling convention and a different
// register file. This probe is the test.
// Two passes, because "what does this dimension cost" and "does this dimension
// fit" are different questions and one sweep cannot answer both.
// PASS 1 gives every dimension a generous stack, so each one returns and
// reports what it actually needs. PASS 2 gives each dimension the control
// task's own stack size, so the figures share that denominator and the first
// dimension that does not fit establishes the supported maximum.
// Each dimension runs in ITS OWN task. The high-water mark is a running
// minimum over a task's life, so one task sweeping every dimension would report
// the deepest dimension's figure against all of them. One task per dimension
// makes each figure that dimension's own.
// Dimensions run in INCREASING order and each reports before the next is
// attempted, so a dimension that does not fit leaves a complete record of every
// dimension below it.

#include "stack_probe.h"
#include "stack_probe_dimension.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "freertos/semphr.h"

#include "esp_log.h"
#include "esp_rom_sys.h"

#include <cstdint>
#include <cinttypes>

namespace {

using ctrlpp::probe::kProbeTag;

// Deliberately generous, so that pass one measures what each dimension NEEDS
// rather than only whether it fitted the control task's stack. A dimension that
// overflows here would be reported the same way, by name.
constexpr std::uint32_t kHeadroomStackBytes = 48 * 1024;

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
    ESP_LOGI(kProbeTag, "PASS 1 -- requirement per dimension, %" PRIu32 " byte stack, scalar = float", kHeadroomStackBytes);
    sweep(kHeadroomStackBytes, done);
    ESP_LOGI(kProbeTag, "PASS 1 complete -- every dimension returned");
}

// Expected to end in a named overflow, and that is the point: it establishes
// the supported maximum by reaching the first dimension that does not fit.
void limit_pass(SemaphoreHandle_t done)
{
    ESP_LOGI(kProbeTag, "PASS 2 -- supported maximum at the control task's own %" PRIu32 " byte stack", ctrlpp::kControlTaskStackBytes);
    sweep(ctrlpp::kControlTaskStackBytes, done);
    ESP_LOGI(kProbeTag, "PASS 2 complete -- no dimension overflowed, which was not expected");
}

void stack_probe_task(void *)
{
    SemaphoreHandle_t done = xSemaphoreCreateBinary();
    headroom_pass(done);
    limit_pass(done);
    vTaskDelete(nullptr);
}

}

namespace ctrlpp {

// The driver task itself needs almost no stack: every solve happens in a task
// this one creates, and it only waits on them.
void start_stack_probe()
{
    xTaskCreatePinnedToCore(stack_probe_task, "ctrlpp_probe", 3 * 1024, nullptr, 4, nullptr, 0);
}

}

// Name the task that overflowed instead of letting the board reset unexplained.
// The framework supplies a weak default; this overrides it so a negative result
// is evidence rather than an absence of output.
// A ROM print and not a log line: the log path takes locks and allocates, and a
// task that has just overrun its stack is the last context that should be
// trusted to complete it. `esp_rom_printf` writes the console UART directly and
// is what survives if the log path does not. Its conversions are the ROM
// printer's own, which is why the figure goes through %u.
extern "C" void vApplicationStackOverflowHook(TaskHandle_t, char *name)
{
    esp_rom_printf("STACK OVERFLOW in task '%s' -- the dimension above the last "
                   "reported probe line does not fit %u bytes\n",
                   name != nullptr ? name : "<unnamed>", static_cast<unsigned>(ctrlpp::kControlTaskStackBytes));
}
