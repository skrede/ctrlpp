#include "instruction_counter.h"

namespace ctrlpp::bench::qp_solvers
{

#if defined(__linux__) && defined(__x86_64__)

namespace
{

auto open_counter() -> int
{
    perf_event_attr attr;
    std::memset(&attr, 0, sizeof(attr));
    attr.type = PERF_TYPE_HARDWARE;
    attr.size = sizeof(attr);
    attr.config = PERF_COUNT_HW_INSTRUCTIONS;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    return static_cast<int>(syscall(SYS_perf_event_open, &attr, 0, -1, -1, 0UL));
}

auto read_pmc(std::uint32_t counter) -> std::uint64_t
{
    std::uint32_t low = 0;
    std::uint32_t high = 0;
    __asm__ volatile("rdpmc" : "=a"(low), "=d"(high) : "c"(counter));
    return (static_cast<std::uint64_t>(high) << 32) | low;
}

}

instruction_counter::instruction_counter()
    : fd_(open_counter())
    , page_(nullptr)
{
    if(fd_ == -1)
        return;
    const auto size = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
    void* page = mmap(nullptr, size, PROT_READ, MAP_SHARED, fd_, 0);
    if(page == MAP_FAILED)
        return;
    if(static_cast<const perf_event_mmap_page*>(page)->cap_user_rdpmc == 0)
    {
        munmap(page, size);
        return;
    }
    page_ = page;
}

instruction_counter::~instruction_counter()
{
    if(page_ != nullptr)
        munmap(page_, static_cast<std::size_t>(sysconf(_SC_PAGESIZE)));
    if(fd_ != -1)
        close(fd_);
}

/// The sequence-lock read of a self-monitored counter: the kernel's offset plus the sign-extended
/// hardware count, retried until the page's lock word is unchanged around the read.
auto instruction_counter::read() const -> std::int64_t
{
    if(page_ == nullptr)
        return instructions_unavailable;
    const auto* header = static_cast<volatile const perf_event_mmap_page*>(page_);
    std::uint32_t sequence = 0;
    std::int64_t count = 0;
    do
    {
        sequence = header->lock;
        __asm__ volatile("" ::: "memory");
        const std::uint32_t index = header->index;
        count = header->offset;
        if(index == 0)
            return instructions_unavailable;
        const auto width = static_cast<std::uint32_t>(header->pmc_width);
        auto pmc = static_cast<std::int64_t>(read_pmc(index - 1) << (64 - width));
        count += pmc >> (64 - width);
        __asm__ volatile("" ::: "memory");
    } while(header->lock != sequence);
    return count;
}

#else

instruction_counter::instruction_counter()
    : fd_(-1)
    , page_(nullptr)
{
}

instruction_counter::~instruction_counter() = default;

auto instruction_counter::read() const -> std::int64_t
{
    return instructions_unavailable;
}

#endif

}
