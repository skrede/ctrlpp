#include "telemetry.h"
#include "alloc_sensor.h"

#include "esp_heap_caps.h"

#include <new>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cinttypes>

// The aligned forms are replaced too: an over-aligned new goes through them and
// never through the plain operator, so leaving them out would leave a path the
// counter cannot see.

namespace {

// Integer conversions into a fixed buffer only, because the heap has just
// refused.
void report_refusal(std::size_t size)
{
    char body[ctrlpp::kLineBytes];
    std::snprintf(body, sizeof(body), "heap REFUSED request=%" PRIu32 " largest_free_block=%" PRIu32 " bytes", static_cast<std::uint32_t>(size),
                  static_cast<std::uint32_t>(heap_caps_get_largest_free_block(MALLOC_CAP_DEFAULT)));
    ctrlpp::send_line(ctrlpp::stream_tag::verdict, body);
}

}

void *operator new(std::size_t size)
{
    if(void *pointer = ctrlpp::allocate_for_new(size > 0 ? size : 1))
        return pointer;
    report_refusal(size);
    std::abort();
}

void *operator new[](std::size_t size)
{
    return ::operator new(size);
}

void *operator new(std::size_t size, std::align_val_t alignment)
{
    if(void *pointer = ctrlpp::allocate_aligned_for_new(size > 0 ? size : 1, static_cast<std::size_t>(alignment)))
        return pointer;
    report_refusal(size);
    std::abort();
}

void *operator new[](std::size_t size, std::align_val_t alignment)
{
    return ::operator new(size, alignment);
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

void operator delete(void *pointer, std::align_val_t) noexcept
{
    ctrlpp::release_for_delete(pointer);
}

void operator delete[](void *pointer, std::align_val_t) noexcept
{
    ctrlpp::release_for_delete(pointer);
}

void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept
{
    ctrlpp::release_for_delete(pointer);
}

void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept
{
    ctrlpp::release_for_delete(pointer);
}
