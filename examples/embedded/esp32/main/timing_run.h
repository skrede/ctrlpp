#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TIMING_RUN_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TIMING_RUN_H

#include "cycle_counter.h"

namespace ctrlpp {

// Proves the counter live, then times the instrument's own overhead, the gain
// design and a golden-length block of each family's steps. Each region's first
// execution after reset is its cold figure, taken before anything else runs it,
// so this must run ahead of the control loop and the families' golden runs.
// Interrupts and the scheduler are left as the demonstration runs them.
counter_state run_timing();

// The estimation family's warm block runs for longer than the task watchdog's
// timeout without yielding, and the watchdog's report would then land inside a
// timed sample. Feeding it or yielding inside the block would change what the
// block measures, so the timing image stops the watchdog before it reports and
// the posture line says so; the evidence image keeps the shipped watchdog.
void stop_task_watchdog() noexcept;

}

#endif
