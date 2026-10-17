// The ukf row's no-malloc executable. This unit arms the QR gain decomposition
// on both sides of the boundary that option brings back, and the part units
// tests/unit/CMakeLists.txt generates beside it arm the allocation-free grid in
// estimation_ukf_nomalloc_grid.h, a few points per unit. The harness header
// must stay the FIRST include of every unit.
//
// The negative control proving both mechanisms fire lives in
// estimation_nomalloc_test.cpp, which is a separate executable; the counter is
// process-global, so every executable consuming the harness runs serially.

#include "nomalloc_harness.h"
#include "estimation_ukf_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>


namespace
{

namespace grid = ctrlpp_test::grid;

}


TEST_CASE("ukf with the QR gain decomposition selected allocates at a measurement dimension of 48",
          "[ukf][hardening][nomalloc][grid]")
{
    SECTION("one below the boundary the decomposition brings back")
    {
        auto created = grid::build_ukf_qr<4, 47>();
        REQUIRE(created.has_value());

        grid::require_alloc_free_steady_state(*created, grid::zero_input(),
                                              grid::Vector<double, 47>::Constant(0.1));
    }

    SECTION("at the boundary")
    {
        auto created = grid::build_ukf_qr<4, 48>();
        REQUIRE(created.has_value());

        grid::require_library_allocation_steady_state(*created, grid::zero_input(),
                                                      grid::Vector<double, 48>::Constant(0.1));
    }
}
