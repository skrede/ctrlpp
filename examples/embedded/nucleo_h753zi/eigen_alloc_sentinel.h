#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_EIGEN_ALLOC_SENTINEL_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_EIGEN_ALLOC_SENTINEL_H

// The Eigen half of the board's two allocation sensors. Every translation unit
// that parses Eigen must include this first, ahead of the project's usual
// include order, because the runtime-no-malloc switch and the eigen_assert
// override only take effect if they are defined before Eigen is parsed; a unit
// that misses them also breaks the one-definition rule for Eigen's malloc gate.
// A failed eigen_assert is counted rather than aborting, so the trap survives
// NDEBUG and -fno-exceptions and the run continues to its report.
//
// There is no over-aligned allocation path to intercept: EIGEN_MAX_ALIGN_BYTES
// is pinned to the natural double width, so Eigen never asks this leg for an
// over-aligned block.

#if defined(EIGEN_CORE_H) || defined(EIGEN_CORE_MODULE_H)
    #error "eigen_alloc_sentinel.h must be the first include of the translation unit, before any Eigen or ctrlpp header"
#endif

// The timing image is built from the same sources with the sentinel off, so the
// code it times carries none of these checks.
#if CTRLPP_MCU_ALLOC_SENTINEL

    #define EIGEN_RUNTIME_NO_MALLOC

    #include <atomic>
    #include <cstdint>

namespace ctrlpp::detail {

inline std::atomic<std::uint32_t> eigen_alloc_violations{0};

}

    #define eigen_assert(X)                                                                                                                                                             \
        do                                                                                                                                                                              \
        {                                                                                                                                                                               \
            if(!(X))                                                                                                                                                                    \
                ::ctrlpp::detail::eigen_alloc_violations.fetch_add(1, ::std::memory_order_relaxed);                                                                                     \
        } while(false)

#endif

#endif
