#include "nomalloc_harness.h"

#include "ctrlpp/mpc/detail/osqp_interop.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>

namespace
{

constexpr double inf = std::numeric_limits<double>::infinity();

auto bound_samples() -> Eigen::VectorXd
{
    Eigen::VectorXd v(8);
    v << -inf, inf, -1e40, 1e40, std::numeric_limits<double>::quiet_NaN(), 0.5, -OSQP_INFTY, OSQP_INFTY;
    return v;
}

}

TEST_CASE("the OSQP bound mapping allocates nothing in steady state", "[osqp][nomalloc]")
{
    const Eigen::VectorXd bounds = bound_samples();
    Eigen::VectorXd lower(bounds.size());
    Eigen::VectorXd upper(bounds.size());
    ctrlpp::detail::to_osqp_lower(bounds, lower);
    ctrlpp::detail::to_osqp_upper(bounds, upper);

    std::size_t allocations = 0;
    bool eigen_violation = false;
    {
        ctrlpp_test::scoped_no_malloc guard;
        for(int i = 0; i < 100; ++i)
        {
            ctrlpp::detail::to_osqp_lower(bounds, lower);
            ctrlpp::detail::to_osqp_upper(bounds, upper);
        }
        allocations = guard.allocations();
        eigen_violation = guard.eigen_violation();
    }

    CHECK(allocations == 0);
    CHECK_FALSE(eigen_violation);

    CHECK(lower(0) == -OSQP_INFTY);
    CHECK(lower(1) == inf);
    CHECK(lower(2) == -OSQP_INFTY);
    CHECK(lower(3) == 1e40);
    CHECK(upper(0) == -inf);
    CHECK(upper(1) == OSQP_INFTY);
    CHECK(upper(2) == -1e40);
    CHECK(upper(3) == OSQP_INFTY);
    CHECK(std::isnan(lower(4)));
    CHECK(std::isnan(upper(4)));
    for(Eigen::Index i = 5; i < bounds.size(); ++i)
    {
        CHECK(lower(i) == bounds(i));
        CHECK(upper(i) == bounds(i));
    }
}
