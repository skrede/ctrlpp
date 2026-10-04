#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_SBRK_CEILING_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_SBRK_CEILING_H

#include <cstddef>

namespace ctrlpp {

std::size_t heap_reserve_bytes() noexcept;

std::size_t heap_high_water_bytes() noexcept;

// The break growth the ceiling refused first, or zero if it has refused none.
std::size_t heap_first_refusal_bytes() noexcept;

}

#endif
