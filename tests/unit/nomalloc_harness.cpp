#include "nomalloc_harness.h"

#include <new>
#include <cstdlib>

#if defined(_MSC_VER)
#include <malloc.h>
#endif


namespace ctrlpp_test::detail
{

std::atomic<std::size_t> allocation_count{0};

}

namespace
{

using ctrlpp_test::detail::allocation_count;

void* counted_try_allocate(std::size_t size) noexcept
{
    allocation_count.fetch_add(1, std::memory_order_relaxed);

    // operator new must return a distinct non-null pointer even for a
    // zero-size request, while std::malloc(0) may return a null pointer.
    return std::malloc(size > 0 ? size : 1);
}

void* counted_try_allocate(std::size_t size, std::align_val_t alignment) noexcept
{
    allocation_count.fetch_add(1, std::memory_order_relaxed);

    // std::aligned_alloc requires the size to be a multiple of the alignment.
    const auto align = static_cast<std::size_t>(alignment);
    const std::size_t padded = align * ((size + align - 1) / align);
    const std::size_t request = padded > 0 ? padded : align;
#if defined(_MSC_VER)
    // MSVC's UCRT ships no std::aligned_alloc; _aligned_malloc takes (size,
    // alignment) and its result must be released with _aligned_free.
    return _aligned_malloc(request, align);
#else
    return std::aligned_alloc(align, request);
#endif
}

void* counted_allocate(std::size_t size)
{
    if(void* pointer = counted_try_allocate(size))
        return pointer;
    std::abort();
}

void* counted_allocate(std::size_t size, std::align_val_t alignment)
{
    if(void* pointer = counted_try_allocate(size, alignment))
        return pointer;
    std::abort();
}

// Aligned allocations from counted_allocate must be released with the matching
// deallocator: _aligned_free on MSVC (_aligned_malloc), std::free elsewhere
// (std::aligned_alloc). The non-aligned deletes pair with std::malloc and stay
// on std::free.
void counted_aligned_free(void* pointer) noexcept
{
#if defined(_MSC_VER)
    _aligned_free(pointer);
#else
    std::free(pointer);
#endif
}

}


void* operator new(std::size_t size)
{
    return counted_allocate(size);
}

void* operator new[](std::size_t size)
{
    return counted_allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
    return counted_allocate(size, alignment);
}

void* operator new[](std::size_t size, std::align_val_t alignment)
{
    return counted_allocate(size, alignment);
}

void* operator new(std::size_t size, const std::nothrow_t&) noexcept
{
    return counted_try_allocate(size);
}

void* operator new[](std::size_t size, const std::nothrow_t&) noexcept
{
    return counted_try_allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment, const std::nothrow_t&) noexcept
{
    return counted_try_allocate(size, alignment);
}

void* operator new[](std::size_t size, std::align_val_t alignment, const std::nothrow_t&) noexcept
{
    return counted_try_allocate(size, alignment);
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
    counted_aligned_free(pointer);
}

void operator delete[](void* pointer, std::align_val_t) noexcept
{
    counted_aligned_free(pointer);
}

void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept
{
    counted_aligned_free(pointer);
}

void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept
{
    counted_aligned_free(pointer);
}

void operator delete(void* pointer, const std::nothrow_t&) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, const std::nothrow_t&) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::align_val_t, const std::nothrow_t&) noexcept
{
    counted_aligned_free(pointer);
}

void operator delete[](void* pointer, std::align_val_t, const std::nothrow_t&) noexcept
{
    counted_aligned_free(pointer);
}
