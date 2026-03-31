#include "ctrlpp/dsp/vector_biquad.h"
#include "ctrlpp/types.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

TEST_CASE("vector_biquad process matches scalar biquad element-wise", "[vector_biquad]")
{
    auto vb = vector_biquad<double, 2>::low_pass(10.0, 100.0);
    auto s0 = biquad<double>::low_pass(10.0, 100.0);
    auto s1 = biquad<double>::low_pass(10.0, 100.0);

    for(int step = 0; step < 10; ++step)
    {
        Vector<double, 2> x;
        x << 0.5, -0.3;

        auto y = vb.process(x);
        auto y0 = s0.process(0.5);
        auto y1 = s1.process(-0.3);

        REQUIRE_THAT(y[0], WithinAbs(y0, 1e-12));
        REQUIRE_THAT(y[1], WithinAbs(y1, 1e-12));
    }
}

TEST_CASE("vector_biquad reset zeroes state", "[vector_biquad]")
{
    auto vb = vector_biquad<double, 2>::low_pass(10.0, 100.0);
    auto fresh = vector_biquad<double, 2>::low_pass(10.0, 100.0);

    Vector<double, 2> x;
    x << 1.0, -1.0;

    for(int i = 0; i < 20; ++i)
        vb.process(x);

    vb.reset();

    for(int i = 0; i < 10; ++i)
    {
        auto y1 = vb.process(x);
        auto y2 = fresh.process(x);
        REQUIRE_THAT(y1[0], WithinAbs(y2[0], 1e-12));
        REQUIRE_THAT(y1[1], WithinAbs(y2[1], 1e-12));
    }
}

TEST_CASE("vector_biquad reset to value", "[vector_biquad]")
{
    auto vb = vector_biquad<double, 2>::low_pass(10.0, 100.0);

    Vector<double, 2> val;
    val << 1.0, 2.0;
    vb.reset(val);

    auto y = vb.process(val);
    REQUIRE_THAT(y[0], WithinAbs(1.0, 0.1));
    REQUIRE_THAT(y[1], WithinAbs(2.0, 0.2));
}

TEST_CASE("vector_biquad notch factory", "[vector_biquad]")
{
    auto vb = vector_biquad<double, 2>::notch(25.0, 100.0, 1.0);

    Vector<double, 2> x;
    x << 1.0, 1.0;

    Vector<double, 2> y;
    for(int i = 0; i < 50; ++i)
        y = vb.process(x);

    REQUIRE_THAT(y[0], WithinAbs(1.0, 0.05));
    REQUIRE_THAT(y[1], WithinAbs(1.0, 0.05));
}

TEST_CASE("vector_cascaded_biquad process matches scalar cascaded_biquad", "[vector_biquad]")
{
    auto vcb = make_vector_butterworth<4, 2>(10.0, 100.0);
    auto s0 = make_butterworth<4>(10.0, 100.0);
    auto s1 = make_butterworth<4>(10.0, 100.0);

    for(int step = 0; step < 10; ++step)
    {
        Vector<double, 2> x;
        x << 0.5, -0.3;

        auto y = vcb.process(x);
        auto y0 = s0.process(0.5);
        auto y1 = s1.process(-0.3);

        REQUIRE_THAT(y[0], WithinAbs(y0, 1e-12));
        REQUIRE_THAT(y[1], WithinAbs(y1, 1e-12));
    }
}

TEST_CASE("make_vector_chebyshev1 produces valid filter", "[vector_biquad]")
{
    auto vcb = make_vector_chebyshev1<4, 2>(10.0, 100.0, 1.0);
    auto scalar_cb = make_chebyshev1<4>(10.0, 100.0, 1.0);

    Vector<double, 2> x;
    x << 1.0, 1.0;

    Vector<double, 2> y;
    double y_scalar = 0.0;
    for(int i = 0; i < 500; ++i)
    {
        y = vcb.process(x);
        y_scalar = scalar_cb.process(1.0);
    }

    // DC gain matches scalar Chebyshev1 (not exactly 1.0 due to ripple normalization)
    REQUIRE_THAT(y[0], WithinAbs(y_scalar, 1e-12));
    REQUIRE_THAT(y[1], WithinAbs(y_scalar, 1e-12));
    // Output is in reasonable range for a low-pass filter at DC
    REQUIRE(y[0] > 0.8);
    REQUIRE(y[1] > 0.8);
}

TEST_CASE("vector_discrete_filter concept check", "[vector_biquad]")
{
    static_assert(vector_discrete_filter<vector_biquad<double, 2>, Vector<double, 2>>);
    static_assert(vector_discrete_filter<vector_cascaded_biquad<double, 2, 2>, Vector<double, 2>>);
    SUCCEED();
}
