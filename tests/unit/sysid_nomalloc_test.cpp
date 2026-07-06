// Verify the steady-state recursive system-identification hot paths do zero heap
// allocation on fixed-size templated inputs, using the belt-and-suspenders
// harness from nomalloc_harness.h: a throwing eigen_assert that survives
// -DNDEBUG plus a global allocation counter that catches heap traffic outside
// Eigen's own bookkeeping. The harness header must stay the first include.
//
// Coverage: rls::update and recursive_arx::update. One warm-up update outside
// the armed window flushes lazy one-time instantiation; only the steady-state
// update loop is guarded.

#include "nomalloc_harness.h"

#include "ctrlpp/sysid/rls.h"
#include "ctrlpp/sysid/recursive_arx.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <utility>


namespace
{

template <typename Fn>
std::size_t guarded_allocations(Fn&& fn)
{
    ctrlpp_test::scoped_no_malloc guard;
    std::forward<Fn>(fn)();
    return guard.allocations();
}

}


TEST_CASE("rls update performs zero heap allocation",
          "[sysid][rls][hardening][nomalloc]")
{
    constexpr std::size_t np = 3;
    ctrlpp::rls<double, np> estimator;

    ctrlpp::Vector<double, np> phi;
    phi << 0.3, -0.7, 0.5;
    const double y = 1.23;

    estimator.update(y, phi);

    std::size_t allocations = 0;
    bool ok = true;
    REQUIRE_NOTHROW(allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
        {
            phi(0) = 0.3 + 0.001 * static_cast<double>(i);
            ok = estimator.update(y, phi) && ok;
        }
    }));

    REQUIRE(allocations == 0);
    REQUIRE(ok);
    REQUIRE(estimator.parameters().allFinite());
}

TEST_CASE("recursive_arx update performs zero heap allocation",
          "[sysid][recursive_arx][hardening][nomalloc]")
{
    constexpr std::size_t na = 2;
    constexpr std::size_t nb = 2;
    ctrlpp::recursive_arx<double, na, nb> arx;

    arx.update(0.0, 0.0);

    std::size_t allocations = 0;
    REQUIRE_NOTHROW(allocations = guarded_allocations([&] {
        double y = 0.0;
        double u_prev = 0.0;
        for(int i = 0; i < 256; ++i)
        {
            const double u = (i % 2 == 0) ? 1.0 : -1.0;
            const double y_new = 0.8 * y + 0.5 * u_prev;
            arx.update(y_new, u);
            y = y_new;
            u_prev = u;
        }
    }));

    REQUIRE(allocations == 0);
    REQUIRE(arx.parameters().allFinite());
}
