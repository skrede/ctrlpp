#include "ctrlpp/model/propagate.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <type_traits>

using DSS = ctrlpp::discrete_state_space<double, 2, 1, 1>;
using CSS = ctrlpp::continuous_state_space<double, 2, 1, 1>;

TEST_CASE("discrete_state_space aggregate initialization", "[state_space]")
{
    DSS sys{};
    sys.A << 1.0, 1.0, 0.0, 1.0;
    sys.B << 0.0, 1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    REQUIRE(sys.A(0, 0) == 1.0);
    REQUIRE(sys.A(0, 1) == 1.0);
    REQUIRE(sys.A(1, 0) == 0.0);
    REQUIRE(sys.A(1, 1) == 1.0);
    REQUIRE(sys.B(1, 0) == 1.0);
}

TEST_CASE("continuous_state_space is distinct from discrete_state_space", "[state_space]")
{
    static_assert(!std::is_same_v<CSS, DSS>);
    SUCCEED();
}

TEST_CASE("propagate computes x_next = A*x + B*u", "[state_space]")
{
    // Double integrator: A = [[1,1],[0,1]], B = [[0],[1]]
    DSS sys{};
    sys.A << 1.0, 1.0, 0.0, 1.0;
    sys.B << 0.0, 1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    ctrlpp::Vector<double, 2> x;
    x << 0.0, 0.0;
    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    // Step 1: x_next = A*[0,0] + B*[1] = [0,0] + [0,1] = [0,1]
    auto x1 = ctrlpp::propagate(sys, x, u);
    REQUIRE_THAT(x1[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(x1[1], Catch::Matchers::WithinAbs(1.0, 1e-12));

    // Step 2: x_next = A*[0,1] + B*[1] = [1,1] + [0,1] = [1,2]
    auto x2 = ctrlpp::propagate(sys, x1, u);
    REQUIRE_THAT(x2[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(x2[1], Catch::Matchers::WithinAbs(2.0, 1e-12));
}

TEST_CASE("output computes y = C*x + D*u", "[state_space]")
{
    DSS sys{};
    sys.A << 1.0, 1.0, 0.0, 1.0;
    sys.B << 0.0, 1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    ctrlpp::Vector<double, 2> x;
    x << 3.0, 5.0;
    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    // y = C*x + D*u = [1,0]*[3,5] + [0]*[1] = [3]
    auto y = ctrlpp::output(sys, x, u);
    REQUIRE_THAT(y[0], Catch::Matchers::WithinAbs(3.0, 1e-12));

    // Also works on continuous_state_space
    CSS csys{};
    csys.A << 1.0, 1.0, 0.0, 1.0;
    csys.B << 0.0, 1.0;
    csys.C << 2.0, 0.0;
    csys.D << 1.0;

    // y = C*x + D*u = [2,0]*[3,5] + [1]*[1] = [6] + [1] = [7]
    auto y2 = ctrlpp::output(csys, x, u);
    REQUIRE_THAT(y2[0], Catch::Matchers::WithinAbs(7.0, 1e-12));
}

TEST_CASE("SISO alias template works", "[state_space]")
{
    static_assert(std::is_same_v<
        ctrlpp::siso_discrete_state_space<double, 2>,
        ctrlpp::discrete_state_space<double, 2, 1, 1>>);
    static_assert(std::is_same_v<
        ctrlpp::siso_continuous_state_space<double, 2>,
        ctrlpp::continuous_state_space<double, 2, 1, 1>>);
    SUCCEED();
}
