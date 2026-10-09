#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_ALLOC_SENSOR_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_ALLOC_SENSOR_H

#include <cstddef>
#include <cstdint>

namespace ctrlpp {

// The framework's own free-heap figures across the armed windows, a view
// independent of both sensors: what the windows kept, and the deepest dip below
// a window's starting free size. The dip is summed over every heap and every
// task, so another task's allocation lands in it too.
struct heap_view
{
    std::int32_t retained;
    std::uint32_t dip;
    std::uint32_t windows;
};

struct canary_verdict
{
    bool live;
    std::uint32_t observed;
    std::uint32_t expected;
    std::uint32_t eigen_trips;
    bool heap_live;
    heap_view heap;
    std::uint32_t requested;
    std::uint32_t heap_self;
};

// The view's start call allocates a snapshot of every heap's minimum after the
// window's opening free size is read, and its stop call frees it after the
// closing one, so each window keeps and dips by exactly that block, its
// self-count, when nothing else in the window touched the heap. The self-count
// is measured by an empty window, so it follows the framework rather than a
// figure copied from one version of it.
inline bool heap_view_clean(const heap_view &view, std::uint32_t self) noexcept
{
    const std::uint32_t dip = view.windows == 0 ? 0 : self;
    return std::int64_t{view.retained} == std::int64_t{view.windows} * self && view.dip == dip;
}

// The evidence image carries the sensor. The timing image is built without it,
// keeps arm, disarm and reset as no-ops so the run is the same, and has no
// count to read: its allocation figures are reported as unmeasured, not zero.
inline constexpr bool kAllocSensorPresent = CTRLPP_MCU_ALLOC_SENTINEL != 0;

// The replaced operators' route to the framework's allocator: counted in the
// evidence image, the plain allocator otherwise.
void *allocate_for_new(std::size_t size) noexcept;

void *allocate_aligned_for_new(std::size_t size, std::size_t alignment) noexcept;

void release_for_delete(void *pointer) noexcept;

// Counting is off until armed, so setup allocation is never charged to the
// window being measured. Arming also forbids Eigen allocation, which is what
// makes the Eigen sentinel trip.
void alloc_sensor_arm() noexcept;

void alloc_sensor_disarm() noexcept;

void alloc_sensor_reset() noexcept;

std::uint32_t alloc_sensor_observed() noexcept;

std::uint32_t alloc_sensor_eigen_trips() noexcept;

heap_view alloc_sensor_heap() noexcept;

// The empty window's figure, read when the canary runs.
std::uint32_t alloc_sensor_heap_self() noexcept;

// Measures the heap view's self-count, then makes one deliberate allocation
// through each path the sensor claims to see inside an armed window, and one
// through the C allocator that only the heap view can see. A zero the sensor
// reports later is evidence only if this came back live; the heap view is live
// if the window kept the requested bytes beyond its self-count.
canary_verdict run_alloc_canary() noexcept;

}

#endif
