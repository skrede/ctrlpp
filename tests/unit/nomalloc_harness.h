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
// 2. The global allocation functions (operator new, operator new[], their
//    std::nothrow_t forms, and the matching deletes) are replaced with counting
//    versions. Every form is replaced because a sanitizer runtime supplies its
//    own definition of any form left out, and its memory then reaches the
//    replaced std::free-based deletes as an allocator mismatch. Eigen's internal
//    aligned_malloc calls std::malloc directly and bypasses operator new, so
//    the Eigen-side check covers what the counter cannot see, while the
//    counter covers every non-Eigen heap allocation that Eigen's bookkeeping
//    cannot see. Either mechanism alone can silently false-pass. A hard
//    allocation failure in a throwing form (which may not return null) calls
//    std::abort, since raising an exception is unavailable under
//    -fno-exceptions; the nothrow forms return null as the standard requires.
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
// The replacement allocation functions and the counter they feed are defined
// once, in nomalloc_harness.cpp, which every consuming executable links exactly
// once through the ctrlpp_nomalloc_harness object library. Several translation
// units of one executable may include this header, each as its first include.
// The counter is process-global, so every consuming executable must be
// registered RUN_SERIAL.
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

#include <cstddef>


namespace ctrlpp_test::detail
{

extern std::atomic<std::size_t> allocation_count;

}

namespace ctrlpp_test
{

inline std::size_t alloc_count()
{
    return detail::allocation_count.load(std::memory_order_relaxed);
}

// True once a failed eigen_assert has fired since the last scoped_no_malloc
// construction reset the sentinel. The Eigen-side trap is throw-free, so
// consumers poll this instead of catching an exception.
inline bool eigen_violation()
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

    std::size_t allocations() const
    {
        return alloc_count() - baseline_;
    }

    bool eigen_violation() const
    {
        return ctrlpp_test::eigen_violation();
    }

private:
    std::size_t baseline_;
};

}

#endif
