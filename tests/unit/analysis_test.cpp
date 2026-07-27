#include "ctrlpp/model/analysis.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <algorithm>

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

// ---------------------------------------------------------------------------
// Indeterminate (non-finite) inputs
//
// Each predicate answers whether the system is provably stable, controllable,
// or observable. A matrix carrying a NaN or an infinity proves nothing, so the
// answer is false. Every case below pairs the non-finite input with the finite
// input it was derived from, so the pre-existing answer is pinned alongside the
// new one.
// ---------------------------------------------------------------------------

TEST_CASE("is_stable answers false for a non-finite continuous state matrix")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = -1.0;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;

    CHECK(ctrlpp::is_stable(sys));

    sys.A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(ctrlpp::is_stable(sys));

    sys.A(0, 0) = -std::numeric_limits<double>::infinity();
    CHECK_FALSE(ctrlpp::is_stable(sys));
}

TEST_CASE("is_stable answers false for a non-finite discrete state matrix")
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.5;
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = 0.3;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;

    CHECK(ctrlpp::is_stable(sys));

    sys.A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(ctrlpp::is_stable(sys));

    sys.A(0, 0) = std::numeric_limits<double>::infinity();
    CHECK_FALSE(ctrlpp::is_stable(sys));
}

TEST_CASE("poles reports a non-finite spectrum for a non-finite state matrix")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    sys.A(0, 1) = 0.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -2.0;

    auto p = ctrlpp::poles(sys);

    // Every entry is marked non-finite, so no pole of the returned spectrum can
    // be mistaken for a computed eigenvalue.
    for(const auto& pole : p)
    {
        CHECK(std::isnan(pole.real()));
        CHECK(std::isnan(pole.imag()));
    }
}

TEST_CASE("is_controllable answers false for a non-finite state matrix")
{
    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 2, 1> B;
    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;

    CHECK(ctrlpp::is_controllable<double, 2, 1>(A, B));

    A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(ctrlpp::is_controllable<double, 2, 1>(A, B));
}

TEST_CASE("is_controllable answers false for a non-finite input matrix")
{
    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 2, 1> B;
    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;

    CHECK(ctrlpp::is_controllable<double, 2, 1>(A, B));

    B(0, 0) = std::numeric_limits<double>::infinity();
    CHECK_FALSE(ctrlpp::is_controllable<double, 2, 1>(A, B));
}

TEST_CASE("is_observable answers false for a non-finite state matrix")
{
    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 1, 2> C;
    A << 1.0, 1.0, 0.0, 1.0;
    C << 1.0, 0.0;

    CHECK(ctrlpp::is_observable<double, 2, 1>(A, C));

    A(1, 1) = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(ctrlpp::is_observable<double, 2, 1>(A, C));
}

TEST_CASE("is_observable answers false for a non-finite output matrix")
{
    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 1, 2> C;
    A << 1.0, 1.0, 0.0, 1.0;
    C << 1.0, 0.0;

    CHECK(ctrlpp::is_observable<double, 2, 1>(A, C));

    C(0, 0) = -std::numeric_limits<double>::infinity();
    CHECK_FALSE(ctrlpp::is_observable<double, 2, 1>(A, C));
}

TEST_CASE("is_stable_closed_loop answers false for a non-finite operand")
{
    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 2, 1> B;
    ctrlpp::Matrix<double, 1, 2> K;
    A << 0.9, 0.1, 0.0, 0.8;
    B << 0.0, 1.0;
    K << 0.0, 0.0;

    CHECK(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));

    SECTION("non-finite state matrix")
    {
        A(0, 0) = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));
    }

    SECTION("non-finite input matrix")
    {
        B(1, 0) = std::numeric_limits<double>::infinity();
        CHECK_FALSE(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));
    }

    SECTION("non-finite gain matrix")
    {
        K(0, 1) = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));
    }
}

TEST_CASE("is_stable_closed_loop answers false when the closed-loop matrix overflows")
{
    // Every operand is finite, but the product B*K overflows the representable
    // range, so A - B*K is not finite and nothing is provable from it.
    constexpr double huge = std::numeric_limits<double>::max();

    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 2, 1> B;
    ctrlpp::Matrix<double, 1, 2> K;
    A << 0.9, 0.1, 0.0, 0.8;
    B << huge, huge;
    K << huge, huge;

    REQUIRE(A.allFinite());
    REQUIRE(B.allFinite());
    REQUIRE(K.allFinite());
    CHECK_FALSE(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));
}

TEST_CASE("is_stable_observer answers false for a non-finite operand")
{
    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 2, 1> L;
    ctrlpp::Matrix<double, 1, 2> C;
    A << 0.9, 0.1, 0.0, 0.8;
    L << 0.5, 0.1;
    C << 1.0, 0.0;

    CHECK(ctrlpp::is_stable_observer<double, 2, 1>(A, L, C));

    SECTION("non-finite state matrix")
    {
        A(1, 1) = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(ctrlpp::is_stable_observer<double, 2, 1>(A, L, C));
    }

    SECTION("non-finite observer gain matrix")
    {
        L(0, 0) = -std::numeric_limits<double>::infinity();
        CHECK_FALSE(ctrlpp::is_stable_observer<double, 2, 1>(A, L, C));
    }

    SECTION("non-finite output matrix")
    {
        C(0, 1) = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(ctrlpp::is_stable_observer<double, 2, 1>(A, L, C));
    }
}

TEST_CASE("is_stable_observer answers false when the observer matrix overflows")
{
    // Same overflow argument as the closed-loop case: L*C leaves the
    // representable range even though every operand is inside it.
    constexpr double huge = std::numeric_limits<double>::max();

    ctrlpp::Matrix<double, 2, 2> A;
    ctrlpp::Matrix<double, 2, 1> L;
    ctrlpp::Matrix<double, 1, 2> C;
    A << 0.9, 0.1, 0.0, 0.8;
    L << huge, huge;
    C << huge, huge;

    REQUIRE(A.allFinite());
    REQUIRE(L.allFinite());
    REQUIRE(C.allFinite());
    CHECK_FALSE(ctrlpp::is_stable_observer<double, 2, 1>(A, L, C));
}
