// Verify the steady-state DSP filter hot paths do zero heap allocation on
// fixed-size templated inputs, using the belt-and-suspenders harness from
// nomalloc_harness.h: a throwing eigen_assert that survives -DNDEBUG plus a
// global allocation counter that catches heap traffic outside Eigen's own
// bookkeeping. The harness header must stay the first include of this file.
//
// Coverage: biquad, cascaded_biquad (built via make_butterworth), vector_biquad,
// and fir process(). The factories return ctrlpp::expected and are unwrapped
// outside the armed window; only the steady-state process() loop is guarded.

#include "nomalloc_harness.h"

#include "ctrlpp/dsp/fir.h"
#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/vector_biquad.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <array>
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


TEST_CASE("biquad process performs zero heap allocation",
          "[dsp][biquad][hardening][nomalloc]")
{
    auto filter = ctrlpp::biquad<double>::low_pass(100.0, 1000.0);
    REQUIRE(filter.has_value());
    auto& bq = *filter;

    for(int i = 0; i < 8; ++i)
        bq.process(1.0);

    std::size_t allocations = 0;
    double y = 0.0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            y = bq.process(1.0);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(std::isfinite(y));
}

TEST_CASE("cascaded_biquad process performs zero heap allocation",
          "[dsp][biquad][hardening][nomalloc]")
{
    auto filter = ctrlpp::make_butterworth<4>(100.0, 1000.0);
    REQUIRE(filter.has_value());
    auto& cascade = *filter;

    for(int i = 0; i < 8; ++i)
        cascade.process(1.0);

    std::size_t allocations = 0;
    double y = 0.0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            y = cascade.process(1.0);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(std::isfinite(y));
}

TEST_CASE("vector_biquad process performs zero heap allocation",
          "[dsp][biquad][hardening][nomalloc]")
{
    constexpr std::size_t N = 4;

    auto filter = ctrlpp::vector_biquad<double, N>::low_pass(100.0, 1000.0);
    REQUIRE(filter.has_value());
    auto& vbq = *filter;

    const ctrlpp::Vector<double, N> x = ctrlpp::Vector<double, N>::Constant(1.0);

    for(int i = 0; i < 8; ++i)
        vbq.process(x);

    std::size_t allocations = 0;
    ctrlpp::Vector<double, N> y = ctrlpp::Vector<double, N>::Zero();
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            y = vbq.process(x);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(y.allFinite());
}

TEST_CASE("fir process performs zero heap allocation",
          "[dsp][fir][hardening][nomalloc]")
{
    ctrlpp::fir<double, 5> filter{std::array<double, 5>{0.1, 0.2, 0.4, 0.2, 0.1}};

    for(int i = 0; i < 8; ++i)
        filter.process(1.0);

    std::size_t allocations = 0;
    double y = 0.0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 256; ++i)
            y = filter.process(1.0);
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(std::isfinite(y));
}
