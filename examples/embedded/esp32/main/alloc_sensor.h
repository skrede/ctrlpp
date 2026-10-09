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
};

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

// Makes one deliberate allocation through each path the sensor claims to see
// inside an armed window, and one through the C allocator that only the heap
// view can see. A zero the sensor reports later is evidence only if this came
// back live.
canary_verdict run_alloc_canary() noexcept;

}

#endif
