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
    static_assert(std::is_same_v<ctrlpp::siso_discrete_state_space<double, 2>, ctrlpp::discrete_state_space<double, 2, 1, 1>>);
    static_assert(std::is_same_v<ctrlpp::siso_continuous_state_space<double, 2>, ctrlpp::continuous_state_space<double, 2, 1, 1>>);
    SUCCEED();
}

TEST_CASE("value-initialized state-space models are the zero system", "[state_space]")
{
    static_assert(std::is_aggregate_v<DSS> && std::is_aggregate_v<CSS>);

    const ctrlpp::discrete_state_space<double, 3, 2, 4> dsys{};
    REQUIRE(dsys.A.isZero(0.0));
    REQUIRE(dsys.B.isZero(0.0));
    REQUIRE(dsys.C.isZero(0.0));
    REQUIRE(dsys.D.isZero(0.0));

    const ctrlpp::continuous_state_space<float, 2, 3, 1> csys{};
    REQUIRE(csys.A.isZero(0.0F));
    REQUIRE(csys.B.isZero(0.0F));
    REQUIRE(csys.C.isZero(0.0F));
    REQUIRE(csys.D.isZero(0.0F));
}

TEST_CASE("designated initialization sets the named matrices and zeroes the rest", "[state_space]")
{
    const DSS sys{.A = ctrlpp::Matrix<double, 2, 2>::Identity(), .C = ctrlpp::Matrix<double, 1, 2>::Ones()};
    REQUIRE(sys.A.isIdentity(0.0));
    REQUIRE(sys.B.isZero(0.0));
    REQUIRE((sys.C.array() == 1.0).all());
    REQUIRE(sys.D.isZero(0.0));

    const CSS csys{.B = ctrlpp::Vector<double, 2>::Ones()};
    REQUIRE(csys.A.isZero(0.0));
    REQUIRE((csys.B.array() == 1.0).all());
    REQUIRE(csys.C.isZero(0.0));
    REQUIRE(csys.D.isZero(0.0));
}
