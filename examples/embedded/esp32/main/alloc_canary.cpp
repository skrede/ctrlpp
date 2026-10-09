#include "eigen_alloc_sentinel.h"

#include "alloc_sensor.h"

#include <Eigen/Core>

#include <new>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <algorithm>

#if CTRLPP_MCU_ALLOC_SENTINEL

namespace {

constexpr std::size_t kCanaryBytes     = sizeof(double);
constexpr std::size_t kCanaryAlignment = 2 * __STDCPP_DEFAULT_NEW_ALIGNMENT__;
// A dynamic-size vector always takes heap storage, so one element is enough.
constexpr Eigen::Index kCanaryElements = 1;
// Four blocks of the same size, two of them through the replaced operators.
constexpr std::uint32_t kCanaryBlocks    = 4;
constexpr std::uint32_t kCanaryCounted   = 2;
constexpr std::uint32_t kCanaryRequested = kCanaryBlocks * kCanaryBytes;

std::uint32_t heap_self = 0;

// Hands the address to an opaque consumer so the compiler cannot prove a
// deliberate allocation unused and remove it together with its free.
void escape(void *pointer) noexcept
{
    __asm__ __volatile__("" : : "r"(pointer) : "memory");
}

// An empty window whose dip differs from what it kept saw another task's
// allocation, so its figure is not the self-count.
bool measure_heap_self() noexcept
{
    ctrlpp::alloc_sensor_reset();
    ctrlpp::alloc_sensor_arm();
    ctrlpp::alloc_sensor_disarm();
    const ctrlpp::heap_view empty = ctrlpp::alloc_sensor_heap();
    heap_self                     = static_cast<std::uint32_t>(std::max(empty.retained, std::int32_t{0}));
    return empty.retained >= 0 && heap_self == empty.dip;
}

}

namespace ctrlpp {

std::uint32_t alloc_sensor_heap_self() noexcept
{
    return heap_self;
}

// Exactly one count per allocation through a replaced operator is required:
// fewer means a path is blind, more means one is counted twice. Every block is
// still held when the window closes, so the heap view must have kept at least
// the bytes requested beyond its self-count.
canary_verdict run_alloc_canary() noexcept
{
    const bool isolated = measure_heap_self();
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
    const heap_view heap         = alloc_sensor_heap();
    ::operator delete(aligned_block, std::align_val_t{kCanaryAlignment});
    ::operator delete(cpp_block);
    std::free(c_block);
    const bool heap_live = isolated && heap.retained - static_cast<std::int32_t>(heap_self) >= static_cast<std::int32_t>(kCanaryRequested);
    return {observed == kCanaryCounted && eigen > 0, observed, kCanaryCounted, eigen, heap_live, heap, kCanaryRequested, heap_self};
}

}

#endif
