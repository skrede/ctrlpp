// A strong _sbrk overriding the one --specs=nosys.specs supplies, which advances
// the break with no upper bound, so an over-budget allocation would run on into
// whatever the memory map places next. This one refuses any request past _eheap
// with the C library's failure value.

#include "sbrk_ceiling.h"

#include <cerrno>
#include <cstddef>
#include <cstdint>

// The heap is a byte range delimited by these two linker symbols with no object
// behind it, so its extent and its break can only be reached through their
// addresses.
extern "C" char end[];
extern "C" char _eheap[];

namespace {

std::size_t heap_used        = 0;
std::size_t heap_used_record = 0;
std::size_t first_refusal    = 0;

}

namespace ctrlpp {

std::size_t heap_reserve_bytes() noexcept
{
    return static_cast<std::size_t>(reinterpret_cast<std::uintptr_t>(_eheap) - reinterpret_cast<std::uintptr_t>(end));
}

std::size_t heap_high_water_bytes() noexcept
{
    return heap_used_record;
}

std::size_t heap_first_refusal_bytes() noexcept
{
    return first_refusal;
}

}

extern "C" void *_sbrk(std::ptrdiff_t increment)
{
    const std::size_t previous = heap_used;
    const bool shrinks         = increment < 0;
    // Unsigned negation, so the most negative increment has a magnitude too.
    const std::size_t magnitude = shrinks ? std::size_t{0} - static_cast<std::size_t>(increment) : static_cast<std::size_t>(increment);
    const std::size_t room      = shrinks ? previous : ctrlpp::heap_reserve_bytes() - previous;
    if(magnitude > room)
    {
        if(first_refusal == 0)
            first_refusal = magnitude;
        errno = ENOMEM;
        return reinterpret_cast<void *>(-1);
    }

    heap_used = shrinks ? previous - magnitude : previous + magnitude;
    if(heap_used > heap_used_record)
        heap_used_record = heap_used;
    return end + previous;
}
