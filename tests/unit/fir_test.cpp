#include "ctrlpp/dsp/fir.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

TEST_CASE("fir impulse response matches tap coefficients", "[fir]")
{
    fir f{std::array{0.25, 0.5, 0.25}};

    CHECK_THAT(f.process(1.0), WithinAbs(0.25, 1e-15));
    CHECK_THAT(f.process(0.0), WithinAbs(0.5, 1e-15));
    CHECK_THAT(f.process(0.0), WithinAbs(0.25, 1e-15));
    CHECK_THAT(f.process(0.0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("fir moving average smooths step input", "[fir]")
{
    fir f{std::array{0.2, 0.2, 0.2, 0.2, 0.2}};

    // Step from 0 to 1: outputs ramp up as delay line fills
    CHECK_THAT(f.process(1.0), WithinAbs(0.2, 1e-15));
    CHECK_THAT(f.process(1.0), WithinAbs(0.4, 1e-15));
    CHECK_THAT(f.process(1.0), WithinAbs(0.6, 1e-15));
    CHECK_THAT(f.process(1.0), WithinAbs(0.8, 1e-15));
    CHECK_THAT(f.process(1.0), WithinAbs(1.0, 1e-15));
}

TEST_CASE("fir reset clears delay line", "[fir]")
{
    fir f{std::array{0.25, 0.5, 0.25}};

    // Build up state
    f.process(1.0);
    f.process(1.0);
    f.process(1.0);

    f.reset();

    // Impulse response should match taps again
    CHECK_THAT(f.process(1.0), WithinAbs(0.25, 1e-15));
    CHECK_THAT(f.process(0.0), WithinAbs(0.5, 1e-15));
    CHECK_THAT(f.process(0.0), WithinAbs(0.25, 1e-15));
}

TEST_CASE("fir taps accessor returns original coefficients", "[fir]")
{
    std::array<double, 3> expected{0.25, 0.5, 0.25};
    fir f{expected};

    auto const& t = f.taps();
    CHECK(t[0] == expected[0]);
    CHECK(t[1] == expected[1]);
    CHECK(t[2] == expected[2]);
}

TEST_CASE("fir convolution with known sequence", "[fir]")
{
    // taps={1,-1} computes first difference
    fir f{std::array{1.0, -1.0}};

    // Input: 1, 2, 3, 4
    // Output: 1*1 + (-1)*0 = 1
    //         1*2 + (-1)*1 = 1
    //         1*3 + (-1)*2 = 1
    //         1*4 + (-1)*3 = 1
    CHECK_THAT(f.process(1.0), WithinAbs(1.0, 1e-15));
    CHECK_THAT(f.process(2.0), WithinAbs(1.0, 1e-15));
    CHECK_THAT(f.process(3.0), WithinAbs(1.0, 1e-15));
    CHECK_THAT(f.process(4.0), WithinAbs(1.0, 1e-15));
}

TEST_CASE("fir all-pass single tap", "[fir]")
{
    fir f{std::array{1.0}};

    CHECK_THAT(f.process(3.14), WithinAbs(3.14, 1e-15));
    CHECK_THAT(f.process(-2.7), WithinAbs(-2.7, 1e-15));
    CHECK_THAT(f.process(0.0), WithinAbs(0.0, 1e-15));
}
