#ifndef HPP_GUARD_CTRLPP_TESTS_NOMALLOC_HARNESS_H
#define HPP_GUARD_CTRLPP_TESTS_NOMALLOC_HARNESS_H

// Belt-and-suspenders no-malloc harness for steady-state hot-path tests.
//
// This header MUST be the FIRST include of its translation unit, before any
// Eigen or ctrlpp header, so both defenses are armed before Eigen is parsed:
//
// 1. eigen_assert is redefined to throw std::runtime_error before any Eigen
//    header is included. EIGEN_RUNTIME_NO_MALLOC reports violations through
//    eigen_assert, and the stock assert is elided under -DNDEBUG, so only the
//    throwing form keeps the trap alive in Release builds.
// 2. The global allocation functions (operator new, operator new[], and the
//    matching deletes) are replaced with counting versions. Eigen's internal
//    aligned_malloc calls std::malloc directly and bypasses operator new, so
//    the Eigen-side check covers what the counter cannot see, while the
//    counter covers every non-Eigen heap allocation that Eigen's bookkeeping
//    cannot see. Either mechanism alone can silently false-pass.
//
// Usage contract:
//   1. Construct inputs and run one warm-up call OUTSIDE the guarded window
//      to flush lazy one-time instantiation.
//   2. Construct a ctrlpp_test::scoped_no_malloc guard AFTER the warm-up and
//      run the steady-state calls. Sample guard.allocations() immediately
//      after the guarded calls, before any test framework macro runs inside
//      the window, because framework macros may themselves allocate.
//   3. REQUIRE(guard.allocations() == 0) and REQUIRE_NOTHROW on the guarded
//      calls; a thrown eigen_assert marks an Eigen-side violation.
//
// The allocation counter is process-global: every test consuming this harness
// must be its own translation unit and executable, must include this header
// exactly once, and must be registered RUN_SERIAL.
//
// The include order in this header is load-bearing and intentionally deviates
// from the project ordering: the macros and <stdexcept> must precede the
// Eigen include, and the Eigen include must precede the helpers below.

#if defined(EIGEN_CORE_H) || defined(EIGEN_CORE_MODULE_H)
#error "nomalloc_harness.h must be the first include of the translation unit, before any Eigen or ctrlpp header"
#endif

#define EIGEN_RUNTIME_NO_MALLOC

#include <stdexcept>

#define eigen_assert(X)                                    \
    do                                                     \
    {                                                      \
        if(!(X))                                           \
            throw std::runtime_error("eigen_assert");      \
    } while(false)

#include <Eigen/Dense>

#include <new>
#include <atomic>
#include <cstddef>
#include <cstdlib>


namespace ctrlpp_test
{

namespace detail
{

std::atomic<std::size_t> allocation_count{0};

inline void* counted_allocate(std::size_t size)
{
    allocation_count.fetch_add(1, std::memory_order_relaxed);

    // operator new must return a distinct non-null pointer even for a
    // zero-size request, while std::malloc(0) may return a null pointer.
    if(void* pointer = std::malloc(size > 0 ? size : 1))
        return pointer;
    throw std::bad_alloc{};
}

inline void* counted_allocate(std::size_t size, std::align_val_t alignment)
{
    allocation_count.fetch_add(1, std::memory_order_relaxed);

    // std::aligned_alloc requires the size to be a multiple of the alignment.
    const auto align = static_cast<std::size_t>(alignment);
    const std::size_t padded = align * ((size + align - 1) / align);
    if(void* pointer = std::aligned_alloc(align, padded > 0 ? padded : align))
        return pointer;
    throw std::bad_alloc{};
}

}

}


void* operator new(std::size_t size)
{
    return ctrlpp_test::detail::counted_allocate(size);
}

void* operator new[](std::size_t size)
{
    return ctrlpp_test::detail::counted_allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
    return ctrlpp_test::detail::counted_allocate(size, alignment);
}

void* operator new[](std::size_t size, std::align_val_t alignment)
{
    return ctrlpp_test::detail::counted_allocate(size, alignment);
}

void operator delete(void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::align_val_t) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, std::align_val_t) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept
{
    std::free(pointer);
}


namespace ctrlpp_test
{

[[nodiscard]] inline std::size_t alloc_count()
{
    return detail::allocation_count.load(std::memory_order_relaxed);
}

// RAII no-malloc window: records the allocation-counter baseline and forbids
// Eigen heap allocation on construction, restores Eigen allocation on
// destruction. allocations() reports the heap allocations observed since the
// baseline and is valid both inside the window and after it ends.
class scoped_no_malloc
{
public:
    scoped_no_malloc()
        : baseline_{alloc_count()}
    {
        Eigen::internal::set_is_malloc_allowed(false);
    }

    scoped_no_malloc(const scoped_no_malloc&) = delete;
    scoped_no_malloc& operator=(const scoped_no_malloc&) = delete;

    ~scoped_no_malloc()
    {
        Eigen::internal::set_is_malloc_allowed(true);
    }

    [[nodiscard]] std::size_t allocations() const
    {
        return alloc_count() - baseline_;
    }

private:
    std::size_t baseline_;
};

}

#endif
