#include "ctrlpp/propagate.h"

#include "naive_linalg.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <type_traits>

using DSS = ctrlpp::DiscreteStateSpace<double, 2, 1, 1, NaiveLinalg>;
using CSS = ctrlpp::ContinuousStateSpace<double, 2, 1, 1, NaiveLinalg>;

TEST_CASE("DiscreteStateSpace aggregate initialization", "[state_space]")
{
    DSS sys{
        .A = {{{1.0, 1.0}, {0.0, 1.0}}},
        .B = {{{0.0}, {1.0}}},
        .C = {{{1.0, 0.0}}},
        .D = {{{0.0}}}
    };

    REQUIRE(sys.A[0][0] == 1.0);
    REQUIRE(sys.A[0][1] == 1.0);
    REQUIRE(sys.A[1][0] == 0.0);
    REQUIRE(sys.A[1][1] == 1.0);
    REQUIRE(sys.B[1][0] == 1.0);
}

TEST_CASE("ContinuousStateSpace is distinct from DiscreteStateSpace", "[state_space]")
{
    static_assert(!std::is_same_v<CSS, DSS>);
    SUCCEED();
}

TEST_CASE("propagate computes x_next = A*x + B*u", "[state_space]")
{
    // Double integrator: A = [[1,1],[0,1]], B = [[0],[1]]
    DSS sys{
        .A = {{{1.0, 1.0}, {0.0, 1.0}}},
        .B = {{{0.0}, {1.0}}},
        .C = {{{1.0, 0.0}}},
        .D = {{{0.0}}}
    };

    NaiveLinalg::vector_type<double, 2> x = {0.0, 0.0};
    NaiveLinalg::vector_type<double, 1> u = {1.0};

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
    DSS sys{
        .A = {{{1.0, 1.0}, {0.0, 1.0}}},
        .B = {{{0.0}, {1.0}}},
        .C = {{{1.0, 0.0}}},
        .D = {{{0.0}}}
    };

    NaiveLinalg::vector_type<double, 2> x = {3.0, 5.0};
    NaiveLinalg::vector_type<double, 1> u = {1.0};

    // y = C*x + D*u = [1,0]*[3,5] + [0]*[1] = [3]
    auto y = ctrlpp::output(sys, x, u);
    REQUIRE_THAT(y[0], Catch::Matchers::WithinAbs(3.0, 1e-12));

    // Also works on ContinuousStateSpace
    CSS csys{
        .A = {{{1.0, 1.0}, {0.0, 1.0}}},
        .B = {{{0.0}, {1.0}}},
        .C = {{{2.0, 0.0}}},
        .D = {{{1.0}}}
    };

    // y = C*x + D*u = [2,0]*[3,5] + [1]*[1] = [6] + [1] = [7]
    auto y2 = ctrlpp::output(csys, x, u);
    REQUIRE_THAT(y2[0], Catch::Matchers::WithinAbs(7.0, 1e-12));
}

TEST_CASE("SISO alias template works", "[state_space]")
{
    static_assert(std::is_same_v<
        ctrlpp::SISODiscreteStateSpace<double, 2, NaiveLinalg>,
        ctrlpp::DiscreteStateSpace<double, 2, 1, 1, NaiveLinalg>>);
    static_assert(std::is_same_v<
        ctrlpp::SISOContinuousStateSpace<double, 2, NaiveLinalg>,
        ctrlpp::ContinuousStateSpace<double, 2, 1, 1, NaiveLinalg>>);
    SUCCEED();
}
