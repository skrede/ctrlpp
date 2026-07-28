#include "hardening_helpers.h"
#include "ctrlpp/sysid/recursive_arx.h"
#include "ctrlpp/sysid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <random>

using Catch::Matchers::WithinAbs;

TEST_CASE("Recursive ARX identifies first-order system")
{
    // True system: y(t) = 0.8*y(t-1) + 0.5*u(t-1)
    constexpr std::size_t NA = 1;
    constexpr std::size_t NB = 1;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for(int t = 0; t < 300; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.8 * y + 0.5 * u_prev;
        arx.update(y_new, u);
        y = y_new;
        u_prev = u;
    }

    auto theta = arx.parameters();
    REQUIRE_THAT(theta(0), WithinAbs(0.8, 0.05));
    REQUIRE_THAT(theta(1), WithinAbs(0.5, 0.05));
}

TEST_CASE("Recursive ARX to_state_space returns companion form")
{
    constexpr std::size_t NA = 1;
    constexpr std::size_t NB = 1;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for(int t = 0; t < 300; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.8 * y + 0.5 * u_prev;
        arx.update(y_new, u);
        y = y_new;
        u_prev = u;
    }

    auto ss = arx.to_state_space();
    // For first-order: A = [a1], B = [b1], C = [1], D = [0]
    REQUIRE_THAT(ss.A(0, 0), WithinAbs(0.8, 0.05));
    REQUIRE_THAT(ss.B(0, 0), WithinAbs(0.5, 0.05));
    REQUIRE_THAT(ss.C(0, 0), WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(ss.D(0, 0), WithinAbs(0.0, 1e-15));
}

TEST_CASE("Recursive ARX state-space simulation matches original response")
{
    constexpr std::size_t NA = 1;
    constexpr std::size_t NB = 1;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    constexpr std::size_t N = 300;
    std::array<double, N> y_data{};
    std::array<double, N> u_data{};

    double y = 0.0;
    double u_prev = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        u_data[t] = u_dist(gen);
        double y_new = 0.8 * y + 0.5 * u_prev;
        arx.update(y_new, u_data[t]);
        y_data[t] = y_new;
        y = y_new;
        u_prev = u_data[t];
    }

    auto ss = arx.to_state_space();

    // Simulate state-space model
    Eigen::Matrix<double, 1, 1> x = Eigen::Matrix<double, 1, 1>::Zero();
    double max_error = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << u_data[t];
        auto y_hat = (ss.C * x + ss.D * u_vec).eval();
        x = (ss.A * x + ss.B * u_vec).eval();
        // Skip first few transient samples
        if(t > 10)
        {
            max_error = std::max(max_error, std::abs(y_hat(0, 0) - y_data[t]));
        }
    }
    REQUIRE(max_error < 0.1);
}

TEST_CASE("Recursive ARX delegates parameters and covariance")
{
    constexpr std::size_t NA = 1;
    constexpr std::size_t NB = 1;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    // Initial parameters should be zero
    auto theta = arx.parameters();
    REQUIRE_THAT(theta(0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(theta(1), WithinAbs(0.0, 1e-15));

    // Covariance should be large diagonal initially
    auto P = arx.covariance();
    REQUIRE(P(0, 0) > 100.0);
    REQUIRE(P(1, 1) > 100.0);
}

TEST_CASE("Recursive ARX second-order system identification")
{
    // True system: y(t) = 1.2*y(t-1) - 0.5*y(t-2) + 0.3*u(t-1) + 0.1*u(t-2)
    // (Stable: poles of z^2 - 1.2z + 0.5 are inside unit circle)
    constexpr std::size_t NA = 2;
    constexpr std::size_t NB = 2;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    std::mt19937 gen(99);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y_prev1 = 0.0;
    double y_prev2 = 0.0;
    double u_prev1 = 0.0;
    double u_prev2 = 0.0;

    for(int t = 0; t < 500; ++t)
    {
        double u = u_dist(gen);
        double y_new = 1.2 * y_prev1 - 0.5 * y_prev2 + 0.3 * u_prev1 + 0.1 * u_prev2;
        arx.update(y_new, u);
        y_prev2 = y_prev1;
        y_prev1 = y_new;
        u_prev2 = u_prev1;
        u_prev1 = u;
    }

    auto theta = arx.parameters();
    REQUIRE_THAT(theta(0), WithinAbs(1.2, 0.1));
    REQUIRE_THAT(theta(1), WithinAbs(-0.5, 0.1));
    REQUIRE_THAT(theta(2), WithinAbs(0.3, 0.1));
    REQUIRE_THAT(theta(3), WithinAbs(0.1, 0.1));
}

TEST_CASE("Recursive ARX initial updates before buffer is full do not crash")
{
    constexpr std::size_t NA = 3;
    constexpr std::size_t NB = 2;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    // Just a few updates -- should not crash
    arx.update(1.0, 0.5);
    arx.update(0.8, 0.3);

    // Should have valid (though not converged) parameters
    auto theta = arx.parameters();
    REQUIRE(std::isfinite(theta(0)));
}

TEST_CASE("Recursive ARX with NB > NA realizes all b-coefficients (max(NA,NB) states)")
{
    // True system: y(t) = a1*y(t-1) + b1*u(t-1) + b2*u(t-2), so NA=1 < NB=2.
    // to_state_space() must produce a max(NA,NB)=2 state realization that keeps
    // b2; a NA-state realization would drop it and mispredict the response.
    constexpr double a1 = 0.7;
    constexpr double b1 = 0.5;
    constexpr double b2 = -0.3;

    constexpr std::size_t NA = 1;
    constexpr std::size_t NB = 2;
    auto arx = ctrlpp::test::constructed(ctrlpp::recursive_arx<double, NA, NB>::create());

    std::mt19937 gen(7);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y_prev1 = 0.0;
    double u_prev1 = 0.0;
    double u_prev2 = 0.0;
    for(int t = 0; t < 800; ++t)
    {
        double u = u_dist(gen);
        double y_new = a1 * y_prev1 + b1 * u_prev1 + b2 * u_prev2;
        arx.update(y_new, u);
        y_prev1 = y_new;
        u_prev2 = u_prev1;
        u_prev1 = u;
    }

    auto ss = arx.to_state_space();

    // Realized state dimension is max(NA, NB) = 2, not NA = 1.
    REQUIRE(ss.A.rows() == 2);
    REQUIRE(ss.A.cols() == 2);
    REQUIRE(ss.B.rows() == 2);

    // The identified b2 coefficient survives into the realization and matches
    // the pulse response h(2) = a1*b1 + b2 of the true plant.
    Eigen::Matrix<double, 2, 1> x = Eigen::Matrix<double, 2, 1>::Zero();
    std::array<double, 4> h{};
    for(std::size_t k = 0; k < h.size(); ++k)
    {
        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << (k == 0 ? 1.0 : 0.0);
        h[k] = (ss.C * x + ss.D * u_vec).eval()(0, 0);
        x = (ss.A * x + ss.B * u_vec).eval();
    }

    REQUIRE_THAT(h[1], WithinAbs(b1, 0.05));
    REQUIRE_THAT(h[2], WithinAbs(a1 * b1 + b2, 0.05));
    REQUIRE_THAT(ss.B(1, 0), WithinAbs(b2, 0.05));
}
