#include "alloc_sensor.h"

#include <atomic>
#include <cerrno>
#include <cstddef>

extern "C"
{
void* __libc_malloc(std::size_t);
void __libc_free(void*);
void* __libc_calloc(std::size_t, std::size_t);
void* __libc_realloc(void*, std::size_t);
void* __libc_memalign(std::size_t, std::size_t);
}

namespace
{

constinit std::atomic<std::int64_t> allocations{0};
constinit std::atomic<bool> armed{false};

/// Forwarding to glibc's own entry points avoids the dlsym bootstrap recursion an RTLD_NEXT
/// interposer would hit.
void count_allocation() noexcept
{
    if(armed.load(std::memory_order_relaxed))
        allocations.fetch_add(1, std::memory_order_relaxed);
}

}

extern "C" void* malloc(std::size_t size)
{
    count_allocation();
    return __libc_malloc(size);
}

extern "C" void* calloc(std::size_t count, std::size_t size)
{
    count_allocation();
    return __libc_calloc(count, size);
}

extern "C" void* realloc(void* pointer, std::size_t size)
{
    count_allocation();
    return __libc_realloc(pointer, size);
}

extern "C" void free(void* pointer)
{
    __libc_free(pointer);
}

extern "C" int posix_memalign(void** out, std::size_t alignment, std::size_t size)
{
    count_allocation();
    void* pointer = __libc_memalign(alignment, size);
    if(pointer == nullptr)
        return ENOMEM;
    *out = pointer;
    return 0;
}

extern "C" void* aligned_alloc(std::size_t alignment, std::size_t size)
{
    count_allocation();
    return __libc_memalign(alignment, size);
}

extern "C" void* memalign(std::size_t alignment, std::size_t size)
{
    count_allocation();
    return __libc_memalign(alignment, size);
}

namespace ctrlpp::bench::qp_solvers
{

auto alloc_sensor_present() -> bool
{
    return true;
}

void arm_alloc_sensor()
{
    allocations.store(0, std::memory_order_relaxed);
    armed.store(true, std::memory_order_relaxed);
}

auto disarm_alloc_sensor() -> std::int64_t
{
    armed.store(false, std::memory_order_relaxed);
    return allocations.load(std::memory_order_relaxed);
}

}
