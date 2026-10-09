#include "eigen_alloc_sentinel.h"

#include "alloc_sensor.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_heap_caps.h"

#include <Eigen/Core>

#include <new>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <algorithm>

// The host harness's pairing, not the other board's. This leg's default build
// does not wrap the C allocation family, so the replaced operator new is the
// primary counter and the Eigen sentinel covers the allocations Eigen makes
// through the C allocator, which the counter cannot see.

namespace {

#if CTRLPP_MCU_ALLOC_SENTINEL

// Counted per task, and a null handle is the disarmed state: other tasks,
// preempting this one or running on the other core, allocate on their own
// account, and charging those to the window would misattribute them.
std::atomic<TaskHandle_t> armed_task{nullptr};
std::atomic<std::uint32_t> allocations{0};

std::size_t free_at_arm  = 0;
ctrlpp::heap_view window = {0, 0};

constexpr std::size_t kCanaryBytes     = sizeof(double);
constexpr std::size_t kCanaryAlignment = 2 * __STDCPP_DEFAULT_NEW_ALIGNMENT__;
// A dynamic-size vector always takes heap storage, so one element is enough.
constexpr Eigen::Index kCanaryElements = 1;
// Four blocks of the same size, two of them through the replaced operators.
constexpr std::uint32_t kCanaryBlocks    = 4;
constexpr std::uint32_t kCanaryCounted   = 2;
constexpr std::uint32_t kCanaryRequested = kCanaryBlocks * kCanaryBytes;

void count_allocation() noexcept
{
    const TaskHandle_t task = armed_task.load(std::memory_order_relaxed);
    if(task != nullptr && task == xTaskGetCurrentTaskHandle())
        allocations.fetch_add(1, std::memory_order_relaxed);
}

std::size_t free_bytes() noexcept
{
    return heap_caps_get_free_size(MALLOC_CAP_DEFAULT);
}

// Hands the address to an opaque consumer so the compiler cannot prove a
// deliberate allocation unused and remove it together with its free.
void escape(void *pointer) noexcept
{
    __asm__ __volatile__("" : : "r"(pointer) : "memory");
}

#endif

}

namespace ctrlpp {

#if CTRLPP_MCU_ALLOC_SENTINEL

void *allocate_for_new(std::size_t size) noexcept
{
    count_allocation();
    return std::malloc(size);
}

void *allocate_aligned_for_new(std::size_t size, std::size_t alignment) noexcept
{
    count_allocation();
    return heap_caps_aligned_alloc(alignment, size, MALLOC_CAP_DEFAULT);
}

void alloc_sensor_arm() noexcept
{
    Eigen::internal::set_is_malloc_allowed(false);
    free_at_arm = free_bytes();
    heap_caps_monitor_local_minimum_free_size_start();
    armed_task.store(xTaskGetCurrentTaskHandle(), std::memory_order_relaxed);
}

void alloc_sensor_disarm() noexcept
{
    armed_task.store(nullptr, std::memory_order_relaxed);
    const std::size_t minimum = heap_caps_get_minimum_free_size(MALLOC_CAP_DEFAULT);
    const std::size_t now     = free_bytes();
    heap_caps_monitor_local_minimum_free_size_stop();
    window.retained += static_cast<std::int32_t>(free_at_arm) - static_cast<std::int32_t>(now);
    window.dip = std::max(window.dip, static_cast<std::uint32_t>(free_at_arm - std::min(minimum, free_at_arm)));
    Eigen::internal::set_is_malloc_allowed(true);
}

void alloc_sensor_reset() noexcept
{
    allocations.store(0, std::memory_order_relaxed);
    detail::eigen_alloc_violations.store(0, std::memory_order_relaxed);
    window = {0, 0};
}

std::uint32_t alloc_sensor_observed() noexcept
{
    return allocations.load(std::memory_order_relaxed);
}

std::uint32_t alloc_sensor_eigen_trips() noexcept
{
    return detail::eigen_alloc_violations.load(std::memory_order_relaxed);
}

heap_view alloc_sensor_heap() noexcept
{
    return window;
}

// Exactly one count per allocation through a replaced operator is required:
// fewer means a path is blind, more means one is counted twice. Every block is
// still held when the window closes, so the heap view must have kept at least
// the bytes requested.
canary_verdict run_alloc_canary() noexcept
{
    alloc_sensor_reset();
    alloc_sensor_arm();
    void *c_block               = std::malloc(kCanaryBytes);
    void *cpp_block             = ::operator new(kCanaryBytes);
    void *aligned_block         = ::operator new(kCanaryBytes, std::align_val_t{kCanaryAlignment});
    Eigen::VectorXd eigen_block = Eigen::VectorXd::Ones(kCanaryElements);
    escape(c_block);
    escape(cpp_block);
    escape(aligned_block);
    escape(eigen_block.data());
    alloc_sensor_disarm();

    const std::uint32_t observed = alloc_sensor_observed();
    const std::uint32_t eigen    = alloc_sensor_eigen_trips();
    ::operator delete(aligned_block, std::align_val_t{kCanaryAlignment});
    ::operator delete(cpp_block);
    std::free(c_block);
    const bool heap_live = window.retained >= static_cast<std::int32_t>(kCanaryRequested);
    return {observed == kCanaryCounted && eigen > 0, observed, kCanaryCounted, eigen, heap_live, window, kCanaryRequested};
}

#else

void *allocate_for_new(std::size_t size) noexcept
{
    return std::malloc(size);
}

void *allocate_aligned_for_new(std::size_t size, std::size_t alignment) noexcept
{
    return heap_caps_aligned_alloc(alignment, size, MALLOC_CAP_DEFAULT);
}

void alloc_sensor_arm() noexcept
{
}

void alloc_sensor_disarm() noexcept
{
}

void alloc_sensor_reset() noexcept
{
}

#endif

// The framework's free accepts a block from either allocation call.
void release_for_delete(void *pointer) noexcept
{
    std::free(pointer);
}

}
