#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_INSTRUCTION_COUNTER_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_INSTRUCTION_COUNTER_H

#include <cstdint>

#if defined(__linux__)
    #include <cstring>
    #include <unistd.h>
    #include <sys/mman.h>
    #include <sys/ioctl.h>
    #include <sys/syscall.h>
    #include <linux/perf_event.h>
#endif

namespace ctrlpp::bench::qp_solvers
{

inline constexpr std::int64_t instructions_unavailable = -1;

/// Retired user-space instructions of the calling thread (perf_event_open, PERF_COUNT_HW_INSTRUCTIONS,
/// exclude_kernel), after the RAII counter in argmin's benchmarks/perf_instruction_counter.h. Reads
/// happen between the timed phases of a tick, so they go through the mapped event page with rdpmc, the
/// self-monitoring sequence of the perf_event_open(2) manual page, instead of a read() system call;
/// where rdpmc is unavailable the counter reports itself unavailable rather than put a system call
/// inside a timed span.
class instruction_counter
{
public:
    instruction_counter();

    ~instruction_counter();

    instruction_counter(const instruction_counter&) = delete;
    instruction_counter& operator=(const instruction_counter&) = delete;

    auto available() const -> bool
    {
        return page_ != nullptr;
    }

    auto read() const -> std::int64_t;

private:
    int fd_;
    void* page_;
};

}

#endif
