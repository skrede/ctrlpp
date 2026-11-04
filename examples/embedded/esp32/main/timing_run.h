#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TIMING_RUN_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_TIMING_RUN_H

#include "cycle_counter.h"

namespace ctrlpp
{

// Each region's first execution after reset is its cold figure, so this runs ahead of the control loop and
// the families' golden runs; interrupts and the scheduler are left as the demonstration runs them.
counter_state run_timing();

// The estimation family's warm block outlasts the task watchdog's timeout without yielding, and feeding
// the watchdog or yielding inside the block would change what it measures, so the timing image stops the
// watchdog (the posture line says so) while the evidence image keeps the shipped one.
void stop_task_watchdog() noexcept;

}

#endif
