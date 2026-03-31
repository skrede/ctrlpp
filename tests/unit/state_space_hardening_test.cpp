#include "hardening_helpers.h"

#include "ctrlpp/model/state_space.h"
#include "ctrlpp/types.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

// ── State space hardening ──────────────────────────────────────────────────────

TEST_CASE("State space with all-zero matrices", "[state_space][hardening][negative]")
{
    auto A = ctrlpp::test::zero_matrix<double, 2, 2>();
    auto B = ctrlpp::test::zero_matrix<double, 2, 1>();
    auto C = ctrlpp::test::zero_matrix<double, 1, 2>();
    auto D = ctrlpp::test::zero_matrix<double, 1, 1>();

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{A, B, C, D};

    // Propagate: x(k+1) = A*x(k) + B*u(k) = 0
    ctrlpp::Vector<double, 2> x = ctrlpp::Vector<double, 2>::Ones();
    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    auto x_next = (sys.A * x + sys.B * u).eval();
    REQUIRE_THAT(x_next(0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(x_next(1), WithinAbs(0.0, 1e-15));

    auto y = (sys.C * x + sys.D * u).eval();
    REQUIRE_THAT(y(0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("State space with NaN in system matrices", "[state_space][hardening][negative]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    auto B = ctrlpp::Matrix<double, 2, 1>::Identity();
    auto C = ctrlpp::Matrix<double, 1, 2>::Identity();
    auto D = ctrlpp::test::zero_matrix<double, 1, 1>();

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{A, B, C, D};

    ctrlpp::Vector<double, 2> x = ctrlpp::Vector<double, 2>::Ones();
    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    auto x_next = (sys.A * x + sys.B * u).eval();
    CHECK(std::isnan(x_next(0)));
    CHECK(std::isnan(x_next(1)));
}

TEST_CASE("Identity system propagation matches analytical", "[state_space][hardening][precision]")
{
    // A=I, B=I, C=I, D=0 -> x(k+1) = x(k) + u(k), y(k) = x(k)
    auto A = ctrlpp::Matrix<double, 2, 2>::Identity();
    auto B = ctrlpp::Matrix<double, 2, 2>::Identity();
    auto C = ctrlpp::Matrix<double, 2, 2>::Identity();
    auto D = ctrlpp::test::zero_matrix<double, 2, 2>();

    ctrlpp::discrete_state_space<double, 2, 2, 2> sys{A, B, C, D};

    ctrlpp::Vector<double, 2> x;
    x << 1.0, 2.0;

    ctrlpp::Vector<double, 2> u;
    u << 0.5, -0.3;

    // Step 1: x_next = x + u = [1.5, 1.7]
    auto x_next = (sys.A * x + sys.B * u).eval();
    REQUIRE_THAT(x_next(0), WithinAbs(1.5, 1e-15));
    REQUIRE_THAT(x_next(1), WithinAbs(1.7, 1e-15));

    // y = x
    auto y = (sys.C * x + sys.D * u).eval();
    REQUIRE_THAT(y(0), WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(y(1), WithinAbs(2.0, 1e-15));

    // Step 2: x_next2 = x_next + u = [2.0, 1.4]
    auto x_next2 = (sys.A * x_next + sys.B * u).eval();
    REQUIRE_THAT(x_next2(0), WithinAbs(2.0, 1e-15));
    REQUIRE_THAT(x_next2(1), WithinAbs(1.4, 1e-15));
}

TEST_CASE("Near-singular system from hardening_helpers", "[state_space][hardening][robustness]")
{
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(1e-10);

    ctrlpp::Vector<double, 2> x;
    x << 1.0, 1.0;

    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    // Propagate several steps
    for (int i = 0; i < 100; ++i) {
        x = (sys.A * x + sys.B * u).eval();
    }

    // Near-singular A means first state decays to near-zero, second stays stable
    REQUIRE(std::isfinite(x(0)));
    REQUIRE(std::isfinite(x(1)));

    // First state should have decayed (A(0,0) = 1e-10)
    REQUIRE(std::abs(x(0)) < 2.0); // bounded, not growing
}
