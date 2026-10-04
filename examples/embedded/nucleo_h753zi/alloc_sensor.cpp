#include "eigen_alloc_sentinel.h"

#include "alloc_sensor.h"
#include "sbrk_ceiling.h"

#include <Eigen/Core>

#include <new>
#include <atomic>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cinttypes>

// -Wl,--wrap routes every call that names the C allocation family, from any
// object in the link, to the __wrap_ definitions below and exposes the C
// library's originals as __real_. Calls the C library makes through its
// reentrant entry points (_malloc_r, _calloc_r: the stdio buffer and the float
// formatter's big integers) never name these symbols and are not seen.
extern "C" void *__real_malloc(std::size_t size);
extern "C" void *__real_calloc(std::size_t count, std::size_t size);
extern "C" void *__real_realloc(void *pointer, std::size_t size);
extern "C" void __real_free(void *pointer);

namespace {

std::atomic<bool> armed{false};
std::atomic<std::uint32_t> allocations{0};

constexpr std::size_t kCanaryBytes = sizeof(double);
// A dynamic-size vector always takes heap storage, so one element is enough.
constexpr Eigen::Index kCanaryElements     = 1;
constexpr std::uint32_t kCanaryAllocations = 3;

void count_allocation() noexcept
{
    if(armed.load(std::memory_order_relaxed))
        allocations.fetch_add(1, std::memory_order_relaxed);
}

// Hands the address to an opaque consumer so the compiler cannot prove a
// deliberate allocation unused and remove it together with its free.
void escape(void *pointer) noexcept
{
    __asm__ __volatile__("" : : "r"(pointer) : "memory");
}

}

extern "C" void *__wrap_malloc(std::size_t size)
{
    count_allocation();
    return __real_malloc(size);
}

extern "C" void *__wrap_calloc(std::size_t count, std::size_t size)
{
    count_allocation();
    return __real_calloc(count, size);
}

extern "C" void *__wrap_realloc(void *pointer, std::size_t size)
{
    count_allocation();
    return __real_realloc(pointer, size);
}

extern "C" void __wrap_free(void *pointer)
{
    __real_free(pointer);
}

// Forwarded to __real_malloc and never to malloc: inside this image malloc IS
// the wrapped, counting one, so going through it would count every C++
// allocation twice.
void *operator new(std::size_t size)
{
    count_allocation();
    if(void *pointer = __real_malloc(size > 0 ? size : 1))
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
    __real_free(pointer);
}

void operator delete[](void *pointer) noexcept
{
    __real_free(pointer);
}

void operator delete(void *pointer, std::size_t) noexcept
{
    __real_free(pointer);
}

void operator delete[](void *pointer, std::size_t) noexcept
{
    __real_free(pointer);
}

namespace ctrlpp {

void alloc_sensor_arm() noexcept
{
    Eigen::internal::set_is_malloc_allowed(false);
    armed.store(true, std::memory_order_relaxed);
}

void alloc_sensor_disarm() noexcept
{
    armed.store(false, std::memory_order_relaxed);
    Eigen::internal::set_is_malloc_allowed(true);
}

void alloc_sensor_reset() noexcept
{
    allocations.store(0, std::memory_order_relaxed);
    detail::eigen_alloc_violations.store(0, std::memory_order_relaxed);
}

std::uint32_t alloc_sensor_observed() noexcept
{
    return allocations.load(std::memory_order_relaxed);
}

std::uint32_t alloc_sensor_eigen_trips() noexcept
{
    return detail::eigen_alloc_violations.load(std::memory_order_relaxed);
}

// Exactly one count per deliberate allocation is required: fewer means a path
// is blind, more means one is counted twice.
canary_verdict run_alloc_canary() noexcept
{
    alloc_sensor_reset();
    alloc_sensor_arm();
    void *c_block   = std::malloc(kCanaryBytes);
    void *cpp_block = ::operator new(kCanaryBytes);
    escape(c_block);
    escape(cpp_block);
    {
        Eigen::VectorXd eigen_block = Eigen::VectorXd::Ones(kCanaryElements);
        escape(eigen_block.data());
    }
    alloc_sensor_disarm();

    const std::uint32_t observed = alloc_sensor_observed();
    const std::uint32_t eigen    = alloc_sensor_eigen_trips();
    ::operator delete(cpp_block);
    std::free(c_block);
    return {observed == kCanaryAllocations && eigen > 0, observed, kCanaryAllocations, eigen};
}

}
