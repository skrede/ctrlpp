#ifndef HPP_GUARD_CTRLPP_TESTS_NOMALLOC_HARNESS_H
#define HPP_GUARD_CTRLPP_TESTS_NOMALLOC_HARNESS_H

// Belt-and-suspenders no-malloc harness for steady-state hot-path tests.
//
// This header MUST be the FIRST include of its translation unit, before any
// Eigen or ctrlpp header, so both defenses are armed before Eigen is parsed:
//
// 1. eigen_assert is redefined to set a pollable sentinel
//    (ctrlpp_test::detail::eigen_alloc_violation) instead of raising an
//    exception, so the trap survives -DNDEBUG and also compiles under
//    -fno-exceptions -fno-rtti. EIGEN_RUNTIME_NO_MALLOC reports violations
//    through eigen_assert; a failed predicate stores true into the sentinel
//    and lets execution continue, while the allocation counter still records
//    the allocation, so both defenses stay live.
// 2. The global allocation functions (operator new, operator new[], and the
//    matching deletes) are replaced with counting versions. Eigen's internal
//    aligned_malloc calls std::malloc directly and bypasses operator new, so
//    the Eigen-side check covers what the counter cannot see, while the
//    counter covers every non-Eigen heap allocation that Eigen's bookkeeping
//    cannot see. Either mechanism alone can silently false-pass. A hard
//    allocation failure (operator new may not return null) calls std::abort:
//    raising an exception is unavailable under -fno-exceptions.
//
// Usage contract:
//   1. Construct inputs and run one warm-up call OUTSIDE the guarded window
//      to flush lazy one-time instantiation.
//   2. Construct a ctrlpp_test::scoped_no_malloc guard AFTER the warm-up and
//      run the steady-state calls. Sample guard.allocations() immediately
//      after the guarded calls, before any test framework macro runs inside
//      the window, because framework macros may themselves allocate.
//   3. After the window, assert guard.allocations() == 0 AND
//      !guard.eigen_violation() (equivalently ctrlpp_test::eigen_violation());
//      a set sentinel marks an Eigen-side violation.
//
// The allocation counter is process-global: every test consuming this harness
// must be its own translation unit and executable, must include this header
// exactly once, and must be registered RUN_SERIAL.
//
// The include order in this header is load-bearing and intentionally deviates
// from the project ordering: the sentinel and the eigen_assert macro must
// precede the Eigen include, and the Eigen include must precede the helpers
// below.

#if defined(EIGEN_CORE_H) || defined(EIGEN_CORE_MODULE_H)
#error "nomalloc_harness.h must be the first include of the translation unit, before any Eigen or ctrlpp header"
#endif

#define EIGEN_RUNTIME_NO_MALLOC

#include <atomic>

namespace ctrlpp_test::detail
{

inline std::atomic<bool> eigen_alloc_violation{false};

}

#define eigen_assert(X)                                            \
    do                                                             \
    {                                                              \
        if(!(X))                                                   \
            ::ctrlpp_test::detail::eigen_alloc_violation.store(    \
                true, ::std::memory_order_relaxed);                \
    } while(false)

#include <Eigen/Dense>

#include <new>
#include <cstddef>
#include <cstdlib>

#if defined(_MSC_VER)
#include <malloc.h>
#endif


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
    std::abort();
}

inline void* counted_allocate(std::size_t size, std::align_val_t alignment)
{
    allocation_count.fetch_add(1, std::memory_order_relaxed);

    // std::aligned_alloc requires the size to be a multiple of the alignment.
    const auto align = static_cast<std::size_t>(alignment);
    const std::size_t padded = align * ((size + align - 1) / align);
    const std::size_t request = padded > 0 ? padded : align;
#if defined(_MSC_VER)
    // MSVC's UCRT ships no std::aligned_alloc; _aligned_malloc takes (size,
    // alignment) and its result must be released with _aligned_free.
    if(void* pointer = _aligned_malloc(request, align))
        return pointer;
#else
    if(void* pointer = std::aligned_alloc(align, request))
        return pointer;
#endif
    std::abort();
}

// Aligned allocations from counted_allocate must be released with the matching
// deallocator: _aligned_free on MSVC (_aligned_malloc), std::free elsewhere
// (std::aligned_alloc). The non-aligned deletes pair with std::malloc and stay
// on std::free.
inline void counted_aligned_free(void* pointer) noexcept
{
#if defined(_MSC_VER)
    _aligned_free(pointer);
#else
    std::free(pointer);
#endif
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
    ctrlpp_test::detail::counted_aligned_free(pointer);
}

void operator delete[](void* pointer, std::align_val_t) noexcept
{
    ctrlpp_test::detail::counted_aligned_free(pointer);
}

void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept
{
    ctrlpp_test::detail::counted_aligned_free(pointer);
}

void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept
{
    ctrlpp_test::detail::counted_aligned_free(pointer);
}


namespace ctrlpp_test
{

[[nodiscard]] inline std::size_t alloc_count()
{
    return detail::allocation_count.load(std::memory_order_relaxed);
}

// True once a failed eigen_assert has fired since the last scoped_no_malloc
// construction reset the sentinel. The Eigen-side trap is throw-free, so
// consumers poll this instead of catching an exception.
[[nodiscard]] inline bool eigen_violation()
{
    return detail::eigen_alloc_violation.load(std::memory_order_relaxed);
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
        detail::eigen_alloc_violation.store(false, std::memory_order_relaxed);
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

    [[nodiscard]] bool eigen_violation() const
    {
        return ctrlpp_test::eigen_violation();
    }

private:
    std::size_t baseline_;
};

}

#endif
