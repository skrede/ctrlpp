#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_POSTURE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_POSTURE_H

#include "reset_reason.h"

namespace ctrlpp {

// Two lines: the build and clock the figures were taken under, and the
// scheduling they ran beside. Read from the calling task, which is the one the
// figures time.
void emit_posture(const reset_cause &cause);

}

#endif
