#ifndef HPP_GUARD_CTRLPP_TESTS_ESTIMATION_EKF_NOMALLOC_GRID_H
#define HPP_GUARD_CTRLPP_TESTS_ESTIMATION_EKF_NOMALLOC_GRID_H

// The ekf row's no-allocation dimension grid: every configuration the real-time
// safety matrix prints a per-configuration figure for is armed with both
// mechanisms live, so the heap claim and the published stack figures cover the
// same surface. estimation_nomalloc_grid.h states the boundary the grid found,
// where it comes from, and what the grid does and does not cover.

#include "estimation_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>


namespace ctrlpp_test::ekf_grid
{

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

struct arm_allocation_free
{
    template <std::size_t NX, std::size_t NY>
    void operator()(grid::point<NX, NY>) const
    {
        auto created = grid::build_ekf<NX, NY>();
        REQUIRE(created.has_value());

        grid::require_alloc_free_steady_state(*created, grid::zero_input(),
                                              grid::Vector<double, NY>::Constant(0.1));
    }
};

}

#endif
