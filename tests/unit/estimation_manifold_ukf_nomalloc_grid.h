#ifndef HPP_GUARD_CTRLPP_TESTS_ESTIMATION_MANIFOLD_UKF_NOMALLOC_GRID_H
#define HPP_GUARD_CTRLPP_TESTS_ESTIMATION_MANIFOLD_UKF_NOMALLOC_GRID_H

// The manifold_ukf row's no-allocation dimension grid: every configuration the
// real-time safety matrix prints a per-configuration figure for is armed with
// both mechanisms live, so the heap claim and the published stack figures cover
// the same surface. estimation_nomalloc_grid.h states the boundary the grid
// found, where it comes from, and what the grid does and does not cover.
//
// ONE CALLER AXIS, NOT TWO. This row's state is a rotation, so its state
// dimension is three by construction and is not a caller's to choose. The
// second slot of every point below carries that structural three rather than a
// swept dimension, and a grid giving this row two axes would measure something
// that does not exist. That three is also the width of the gain solve's
// right-hand side, so this row has no single-column exemption from the boundary
// below: every configuration it has is on the blocked path once the measurement
// dimension reaches 48.

#include "estimation_nomalloc_grid.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Geometry>

#include <cstddef>


namespace ctrlpp_test::manifold_ukf_grid
{

// Five measurement rungs (5); five values named only by the supported-maximum
// table -- 1, 12, 20, 28 and 36 (5); and the measurement dimension one below
// the blocked-Householder boundary (1).
constexpr std::size_t allocation_free_points = 5 + 5 + 1;

// The boundary itself and the one instantiation limit this row's single axis
// has.
constexpr std::size_t library_allocating_points = 2;

constexpr std::size_t published_points = allocation_free_points + library_allocating_points;

using allocation_free_grid = grid::over<
    grid::point<2, 3>, grid::point<4, 3>, grid::point<8, 3>, grid::point<16, 3>, grid::point<32, 3>,
    grid::point<1, 3>, grid::point<12, 3>, grid::point<20, 3>, grid::point<28, 3>, grid::point<36, 3>,
    grid::point<47, 3>>;

using library_allocating_grid = grid::over<
    grid::point<48, 3>, grid::point<128, 3>>;

static_assert(allocation_free_grid::size == allocation_free_points);
static_assert(library_allocating_grid::size == library_allocating_points);
static_assert(allocation_free_grid::size + library_allocating_grid::size == published_points,
              "the armed grid must hold exactly the configurations the matrix publishes a figure for");

struct arm_allocation_free
{
    template <std::size_t NY, std::size_t NS>
    void operator()(grid::point<NY, NS>) const
    {
        static_assert(NS == 3, "this row's state is a rotation and its dimension is not a caller's to choose");

        auto created = grid::build_manifold_ukf<NY>();
        REQUIRE(created.has_value());

        grid::require_alloc_free_steady_state(
            *created, grid::body_rate(), grid::gravity_direction<NY>(Eigen::Quaternion<double>::Identity()));
    }
};

}

#endif
