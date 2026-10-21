#ifndef HPP_GUARD_CTRLPP_TESTS_ESTIMATION_UKF_NOMALLOC_GRID_H
#define HPP_GUARD_CTRLPP_TESTS_ESTIMATION_UKF_NOMALLOC_GRID_H

// The ukf row's no-allocation dimension grid: every configuration the real-time
// safety matrix prints a per-configuration figure for is armed with both
// mechanisms live, so the heap claim and the published stack figures cover the
// same surface. estimation_nomalloc_grid.h states the boundary the grid found,
// where it comes from, and what the grid does and does not cover.
//
// THIS ROW IS EXEMPT FROM THAT BOUNDARY BY ITS DEFAULT, AND ONLY BY ITS DEFAULT.
// Its gain decomposition defaults to LDLT, which builds no Householder sequence
// at all, so the whole grid including the corner at 128 states and 128 outputs
// is allocation-free. Selecting the QR option -- a public configuration field --
// puts it on exactly the same blocked path as the rows beside it, and the QR
// case in estimation_ukf_nomalloc_test.cpp walks both sides of that boundary
// rather than leaving the option's behavior asserted from the other rows.

#include "estimation_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>


namespace ctrlpp_test::ukf_grid
{

// Five equal-dimension rungs plus four more on each of the two held-dimension
// lines, the rung where the three lines cross counted once (13); nine values
// named only by the supported-maximum table -- 1 and 12 on the equal line, 1,
// 12, 20 and 28 on the measurement axis, 1, 12 and 20 on the state axis (9); the
// instantiation limit on each axis with the corner where both are largest (3);
// and both sides of the boundary the other rows have, at 47 and at 48, which
// this row's default decomposition does not reach (2). The armed grid holds
// these 27 and one extra below, so it contains them by construction.
constexpr std::size_t published_points = 13 + 9 + 3 + 2;

template <typename... Extras>
using published_plus = grid::over<
    grid::point<2, 2>, grid::point<4, 4>, grid::point<8, 8>, grid::point<16, 16>, grid::point<32, 32>,
    grid::point<4, 2>, grid::point<4, 8>, grid::point<4, 16>, grid::point<4, 32>,
    grid::point<2, 4>, grid::point<8, 4>, grid::point<16, 4>, grid::point<32, 4>,
    grid::point<1, 1>, grid::point<12, 12>,
    grid::point<4, 1>, grid::point<4, 12>, grid::point<4, 20>, grid::point<4, 28>,
    grid::point<1, 4>, grid::point<12, 4>, grid::point<20, 4>,
    grid::point<128, 4>, grid::point<4, 128>, grid::point<128, 128>,
    grid::point<4, 47>, grid::point<4, 48>, Extras...>;

static_assert(published_plus<>::size == published_points,
              "the published list must hold exactly the configurations the matrix publishes a figure for");

// 24 states at 4 outputs is armed although no table prints it: Eigen's
// fixed-size code changes with the dimension, so an armed neighbor does not cover it.
using allocation_free_grid = published_plus<grid::point<24, 4>>;

static_assert(allocation_free_grid::size == published_points + 1);

struct arm_allocation_free
{
    template <std::size_t NX, std::size_t NY>
    void operator()(grid::point<NX, NY>) const
    {
        auto created = grid::build_ukf<NX, NY>();
        REQUIRE(created.has_value());

        grid::require_alloc_free_steady_state(*created, grid::zero_input(),
                                              grid::Vector<double, NY>::Constant(0.1));
    }
};

}

#endif
