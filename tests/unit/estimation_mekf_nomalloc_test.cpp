// The mekf row's no-allocation dimension grid: every configuration the
// real-time safety matrix prints a per-configuration figure for is armed with
// both mechanisms live, so the heap claim and the published stack figures cover
// the same surface. The harness header must stay the FIRST include of this file,
// and estimation_nomalloc_grid.h states the boundary the grid found, where it
// comes from, and what the grid does and does not cover.
//
// THE SELECTOR IS THE BIAS DIMENSION. The error state is the derived
// NE = 3 + NB, which is why this row's instantiation limit on that axis is 125
// where the rows beside it reach 128, and why two cells of its published
// measurement table do not instantiate at all: NB = 2 is refused, so neither
// (2, 2) nor (2, 4) exists to arm. The error state is also the width of the gain
// solve's right-hand side and is never one, so this row has no single-column
// exemption from the boundary below; and that boundary is on the MEASUREMENT
// dimension rather than on the error state, which the point at 125 bias states
// and 4 outputs shows by staying allocation-free.
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

// Eleven of the thirteen configurations its measurement table prints, the other
// two being the NB = 2 cells that do not instantiate (11); four values named
// only by the supported-maximum table -- 12 on the equal line, 20 and 28 on the
// measurement axis, 12 on the bias axis (4); the instantiation limit on the bias
// axis (1); and the measurement dimension one below the blocked-Householder
// boundary (1).
constexpr std::size_t allocation_free_points = 11 + 4 + 1 + 1;

// The boundary itself, the instantiation limit on the measurement axis, and the
// corner where both axes are largest.
constexpr std::size_t library_allocating_points = 3;

constexpr std::size_t published_points = allocation_free_points + library_allocating_points;

using allocation_free_grid = grid::over<
    grid::point<4, 4>, grid::point<8, 8>, grid::point<16, 16>, grid::point<32, 32>,
    grid::point<4, 2>, grid::point<4, 8>, grid::point<4, 16>, grid::point<4, 32>,
    grid::point<8, 4>, grid::point<16, 4>, grid::point<32, 4>,
    grid::point<12, 12>, grid::point<4, 20>, grid::point<4, 28>, grid::point<12, 4>,
    grid::point<125, 4>,
    grid::point<4, 47>>;

using library_allocating_grid = grid::over<
    grid::point<4, 48>, grid::point<4, 128>, grid::point<124, 128>>;

static_assert(allocation_free_grid::size == allocation_free_points);
static_assert(library_allocating_grid::size == library_allocating_points);
static_assert(allocation_free_grid::size + library_allocating_grid::size == published_points,
              "the armed grid must hold exactly the configurations the matrix publishes a figure for");

}


TEST_CASE("mekf predict/update performs zero heap allocation across the published grid",
          "[mekf][hardening][nomalloc][grid]")
{
    const std::size_t armed = allocation_free_grid::arm([]<std::size_t NB, std::size_t NY>(grid::point<NB, NY>) {
        auto created = grid::build_mekf<NB, NY>();
        REQUIRE(created.has_value());

        grid::require_alloc_free_steady_state(
            *created, grid::body_rate(), grid::gravity_direction<NY>(Eigen::Quaternion<double>::Identity()));
    });

    REQUIRE(armed == allocation_free_points);
}

TEST_CASE("mekf allocates in the gain solve at and above a measurement dimension of 48",
          "[mekf][hardening][nomalloc][grid]")
{
    const std::size_t armed = library_allocating_grid::arm([]<std::size_t NB, std::size_t NY>(grid::point<NB, NY>) {
        auto created = grid::build_mekf<NB, NY>();
        REQUIRE(created.has_value());

        grid::require_library_allocation_steady_state(
            *created, grid::body_rate(), grid::gravity_direction<NY>(Eigen::Quaternion<double>::Identity()));
    });

    REQUIRE(armed == library_allocating_points);
}
