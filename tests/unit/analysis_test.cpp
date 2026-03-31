#include "ctrlpp/model/analysis.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

TEST_CASE("poles of stable discrete system")
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.5;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = 0.3;
    sys.B(0, 0) = 1.0;
    sys.B(1, 0) = 0.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    auto p = ctrlpp::poles(sys);

    // Eigenvalues of diagonal matrix are the diagonal elements
    std::array<double, 2> magnitudes{std::abs(p[0]), std::abs(p[1])};
    std::sort(magnitudes.begin(), magnitudes.end());

    CHECK_THAT(magnitudes[0], Catch::Matchers::WithinAbs(0.3, 1e-10));
    CHECK_THAT(magnitudes[1], Catch::Matchers::WithinAbs(0.5, 1e-10));
}

TEST_CASE("is_stable for stable discrete system")
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.5;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = 0.3;
    sys.B(0, 0) = 1.0;
    sys.B(1, 0) = 0.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    CHECK(ctrlpp::is_stable(sys));
}

TEST_CASE("is_stable for unstable discrete system")
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 1.5;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = 0.3;
    sys.B(0, 0) = 1.0;
    sys.B(1, 0) = 0.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    CHECK_FALSE(ctrlpp::is_stable(sys));
}

TEST_CASE("poles of continuous system")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = -1.0;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.B(1, 0) = 0.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    auto p = ctrlpp::poles(sys);

    std::array<double, 2> reals{p[0].real(), p[1].real()};
    std::sort(reals.begin(), reals.end());

    CHECK_THAT(reals[0], Catch::Matchers::WithinAbs(-2.0, 1e-10));
    CHECK_THAT(reals[1], Catch::Matchers::WithinAbs(-1.0, 1e-10));
}

TEST_CASE("is_stable for stable continuous system")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = -1.0;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.B(1, 0) = 0.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    CHECK(ctrlpp::is_stable(sys));
}

TEST_CASE("is_stable for unstable continuous system")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 1.0;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.B(1, 0) = 0.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    CHECK_FALSE(ctrlpp::is_stable(sys));
}
