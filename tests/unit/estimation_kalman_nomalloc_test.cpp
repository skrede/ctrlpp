// The kalman_filter row's no-malloc executable. This unit arms the side of the
// heap claim's boundary where the library allocates, and the part units
// tests/unit/CMakeLists.txt generates beside it arm the allocation-free grid in
// estimation_kalman_nomalloc_grid.h, a few points per unit. The harness header
// must stay the FIRST include of every unit.
//
// The negative control proving both mechanisms fire lives in
// estimation_nomalloc_test.cpp, which is a separate executable; the counter is
// process-global, so every executable consuming the harness runs serially.

#include "nomalloc_harness.h"
#include "estimation_kalman_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>


namespace
{

namespace grid = ctrlpp_test::grid;
namespace row = ctrlpp_test::kalman_grid;

}


TEST_CASE("kalman_filter allocates in the gain solve at and above a measurement dimension of 48",
          "[kalman][hardening][nomalloc][grid]")
{
    const std::size_t armed = row::library_allocating_grid::arm([]<std::size_t NX, std::size_t NY>(grid::point<NX, NY>) {
        auto created = grid::build_kalman<NX, NY>();
        REQUIRE(created.has_value());

        grid::require_library_allocation_steady_state(*created, grid::zero_input(),
                                                      grid::Vector<double, NY>::Constant(0.1));
    });

    REQUIRE(armed == row::library_allocating_points);
}
