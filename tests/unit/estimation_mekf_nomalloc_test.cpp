// The mekf row's no-malloc executable. This unit arms the side of the heap
// claim's boundary where the library allocates, and the part units
// tests/unit/CMakeLists.txt generates beside it arm the allocation-free grid in
// estimation_mekf_nomalloc_grid.h, a few points per unit. The harness header
// must stay the FIRST include of every unit.
//
// The negative control proving both mechanisms fire lives in
// estimation_nomalloc_test.cpp, which is a separate executable; the counter is
// process-global, so every executable consuming the harness runs serially.

#include "nomalloc_harness.h"
#include "estimation_mekf_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Geometry>

#include <cstddef>


namespace
{

namespace grid = ctrlpp_test::grid;
namespace row = ctrlpp_test::mekf_grid;

}


TEST_CASE("mekf allocates in the gain solve at and above a measurement dimension of 48",
          "[mekf][hardening][nomalloc][grid]")
{
    const std::size_t armed = row::library_allocating_grid::arm([]<std::size_t NB, std::size_t NY>(grid::point<NB, NY>) {
        auto created = grid::build_mekf<NB, NY>();
        REQUIRE(created.has_value());

        grid::require_library_allocation_steady_state(
            *created, grid::body_rate(), grid::gravity_direction<NY>(Eigen::Quaternion<double>::Identity()));
    });

    REQUIRE(armed == row::library_allocating_points);
}
