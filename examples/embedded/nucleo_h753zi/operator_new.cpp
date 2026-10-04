#include "alloc_sensor.h"
#include "sbrk_ceiling.h"

#include <new>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cinttypes>

void *operator new(std::size_t size)
{
    if(void *pointer = ctrlpp::allocate_for_new(size > 0 ? size : 1))
        return pointer;
    // Integer conversions only: the float formatter is the one printf path that
    // allocates, and the heap has just refused. The growth is the break's, not
    // this call's size: Eigen reports its own refused malloc by calling here
    // with the largest size there is.
    std::printf("[heap] REFUSED growth=%" PRIu32 " high_water=%" PRIu32 " reserve=%" PRIu32 " bytes\n", static_cast<std::uint32_t>(ctrlpp::heap_first_refusal_bytes()),
                static_cast<std::uint32_t>(ctrlpp::heap_high_water_bytes()), static_cast<std::uint32_t>(ctrlpp::heap_reserve_bytes()));
    std::abort();
}

void *operator new[](std::size_t size)
{
    return ::operator new(size);
}

void operator delete(void *pointer) noexcept
{
    ctrlpp::release_for_delete(pointer);
}

void operator delete[](void *pointer) noexcept
{
    ctrlpp::release_for_delete(pointer);
}

void operator delete(void *pointer, std::size_t) noexcept
{
    ctrlpp::release_for_delete(pointer);
}

void operator delete[](void *pointer, std::size_t) noexcept
{
    ctrlpp::release_for_delete(pointer);
}
