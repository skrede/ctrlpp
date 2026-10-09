#include "eigen_alloc_sentinel.h"

#include "alloc_sensor.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_heap_caps.h"

#include <Eigen/Core>

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
ctrlpp::heap_view window = {0, 0, 0};

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
    ++window.windows;
    Eigen::internal::set_is_malloc_allowed(true);
}

void alloc_sensor_reset() noexcept
{
    allocations.store(0, std::memory_order_relaxed);
    detail::eigen_alloc_violations.store(0, std::memory_order_relaxed);
    window = {0, 0, 0};
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
