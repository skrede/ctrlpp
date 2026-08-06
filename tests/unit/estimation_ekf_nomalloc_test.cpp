// The ekf row's no-allocation dimension grid: every configuration the real-time
// safety matrix prints a per-configuration figure for is armed with both
// mechanisms live, so the heap claim and the published stack figures cover the
// same surface. The harness header must stay the FIRST include of this file, and
// estimation_nomalloc_grid.h states the boundary the grid found, where it comes
// from, and what the grid does and does not cover.
//
// The negative control proving both mechanisms fire lives in
// estimation_nomalloc_test.cpp, which is a separate executable; the counter is
// process-global, so every unit consuming the harness is its own executable and
// is registered to run serially.

#include "nomalloc_harness.h"
#include "estimation_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>


namespace
{

namespace grid = ctrlpp_test::grid;

// Five equal-dimension rungs plus four more on each of the two held-dimension
// lines, the rung where the three lines cross counted once (13); six values
// named only by the supported-maximum table -- 12 on the equal line, 12, 20 and
// 24 on the measurement axis, 20 and 24 on the state axis (6); the instantiation
// limit on the state axis (1); the measurement dimension one below the
// blocked-Householder boundary (1); and that boundary at one state, where the
// right-hand side is a single column and the blocked path is not taken (1).
constexpr std::size_t allocation_free_points = 13 + 6 + 1 + 1 + 1;

// The boundary itself, the instantiation limit on the measurement axis, and the
// corner where both axes are largest. All three carry a multi-column right-hand
// side at a measurement dimension of at least 48.
constexpr std::size_t library_allocating_points = 3;

constexpr std::size_t published_points = allocation_free_points + library_allocating_points;

using allocation_free_grid = grid::over<
    grid::point<2, 2>, grid::point<4, 4>, grid::point<8, 8>, grid::point<16, 16>, grid::point<32, 32>,
    grid::point<4, 2>, grid::point<4, 8>, grid::point<4, 16>, grid::point<4, 32>,
    grid::point<2, 4>, grid::point<8, 4>, grid::point<16, 4>, grid::point<32, 4>,
    grid::point<12, 12>, grid::point<4, 12>, grid::point<4, 20>, grid::point<4, 24>,
    grid::point<20, 4>, grid::point<24, 4>,
    grid::point<128, 4>,
    grid::point<4, 47>,
    grid::point<1, 48>>;

using library_allocating_grid = grid::over<
    grid::point<4, 48>, grid::point<4, 128>, grid::point<128, 128>>;

static_assert(allocation_free_grid::size == allocation_free_points);
static_assert(library_allocating_grid::size == library_allocating_points);
static_assert(allocation_free_grid::size + library_allocating_grid::size == published_points,
              "the armed grid must hold exactly the configurations the matrix publishes a figure for");

}


TEST_CASE("ekf predict/update performs zero heap allocation across the published grid",
          "[ekf][hardening][nomalloc][grid]")
{
    const std::size_t armed = allocation_free_grid::arm([]<std::size_t NX, std::size_t NY>(grid::point<NX, NY>) {
        auto created = grid::build_ekf<NX, NY>();
        REQUIRE(created.has_value());

        grid::require_alloc_free_steady_state(*created, grid::zero_input(),
                                              grid::Vector<double, NY>::Constant(0.1));
    });

    REQUIRE(armed == allocation_free_points);
}

TEST_CASE("ekf allocates in the gain solve at and above a measurement dimension of 48",
          "[ekf][hardening][nomalloc][grid]")
{
    const std::size_t armed = library_allocating_grid::arm([]<std::size_t NX, std::size_t NY>(grid::point<NX, NY>) {
        auto created = grid::build_ekf<NX, NY>();
        REQUIRE(created.has_value());

        grid::require_library_allocation_steady_state(*created, grid::zero_input(),
                                                      grid::Vector<double, NY>::Constant(0.1));
    });

    REQUIRE(armed == library_allocating_points);
}
