#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_RESET_CAUSE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_RESET_CAUSE_H

#include <cstdint>

namespace ctrlpp {

enum class reset_class : std::uint8_t
{
    power_on,
    pin,
    brownout,
    software,
    cpu,
    window_watchdog,
    independent_watchdog,
    d1_standby_exit,
    d2_standby_exit,
    low_power,
    unmatched
};

struct reset_cause
{
    std::uint32_t status;
    reset_class event;
};

// One reset event sets several flags at once, so the flag set is matched whole
// against the events' signatures; a set no single event leaves, such as an
// accumulation that survived a missed clear, is unmatched rather than collapsed.
reset_cause read_reset_cause() noexcept;

// The flags latch and accumulate until removed, so this runs only after the
// cause has been read and printed: the label then describes this run and not
// every event since power-on.
void clear_reset_cause() noexcept;

// Prints the decoded class, every set flag by name and the raw register, as
// fields of the line the caller has open.
void print_reset_cause(const reset_cause &cause);

}

#endif
