#include "ctrlpp/sysid/n4sid.h"
#include "ctrlpp/sysid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <random>

using Catch::Matchers::WithinAbs;

namespace {

// Known 2nd-order discrete state-space system (discretized spring-mass-damper)
// Eigenvalues inside unit circle for stability
constexpr double a11 = 0.9;
constexpr double a12 = 0.1;
constexpr double a21 = -0.1;
constexpr double a22 = 0.85;

struct test_data {
    Eigen::Matrix<double, 1, Eigen::Dynamic> Y;
    Eigen::Matrix<double, 1, Eigen::Dynamic> U;
};

auto generate_data(std::size_t N, double noise_std = 0.0) -> test_data
{
    Eigen::Matrix2d A;
    A << a11, a12, a21, a22;
    Eigen::Vector2d B;
    B << 0.0, 0.1;
    Eigen::RowVector2d C;
    C << 1.0, 0.0;

    Eigen::Matrix<double, 1, Eigen::Dynamic> Y(1, static_cast<Eigen::Index>(N));
    Eigen::Matrix<double, 1, Eigen::Dynamic> U(1, static_cast<Eigen::Index>(N));

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);
    Eigen::Vector2d x = Eigen::Vector2d::Zero();
    for (std::size_t t = 0; t < N; ++t) {
        double u = u_dist(gen);
        U(0, static_cast<Eigen::Index>(t)) = u;

        double y = (C * x)(0);
        if (noise_std > 0.0) {
            std::normal_distribution<double> noise(0.0, noise_std);
            y += noise(gen);
        }
        Y(0, static_cast<Eigen::Index>(t)) = y;

        x = A * x + B * u;
    }

    return {Y, U};
}

}

TEST_CASE("N4SID singular values show clear gap for 2nd-order system") {
    auto [Y, U] = generate_data(1500);

    auto sv = ctrlpp::n4sid_singular_values(Y, U);

    REQUIRE(sv.size() >= 2);
    // First two should be significantly larger than the rest
    REQUIRE(sv(0) > sv(1));
    // Gap: ratio of 2nd to 3rd singular value should be large
    if (sv.size() > 2) {
        double ratio = sv(1) / sv(2);
        REQUIRE(ratio > 5.0);
    }
}

TEST_CASE("N4SID identifies 2nd-order system with good fit") {
    auto [Y, U] = generate_data(1500);

    auto result = ctrlpp::n4sid<2>(Y, U);

    // Simulate identified model on input data and check output match
    REQUIRE(result.metrics.nrmse < 0.1);
}

TEST_CASE("N4SID condition number is finite and positive") {
    auto [Y, U] = generate_data(1500);

    auto result = ctrlpp::n4sid<2>(Y, U);

    REQUIRE(result.condition_number > 0.0);
    REQUIRE(std::isfinite(result.condition_number));
}

TEST_CASE("N4SID VAF > 90 for clean data") {
    auto [Y, U] = generate_data(1500);

    auto result = ctrlpp::n4sid<2>(Y, U);

    REQUIRE(result.metrics.vaf > 90.0);
}

TEST_CASE("N4SID with measurement noise produces degraded but positive VAF") {
    auto [Y, U] = generate_data(1500, 0.05);

    auto result = ctrlpp::n4sid<2>(Y, U);

    REQUIRE(result.metrics.vaf > 0.0);
    // Should still identify something reasonable with mild noise
    REQUIRE(result.metrics.vaf < 100.0);
}
