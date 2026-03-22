#include "ctrlpp/sysid/rls.h"
#include "ctrlpp/sysid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <random>

using Catch::Matchers::WithinAbs;

TEST_CASE("RLS converges to true parameters on known linear system") {
    // y = 2*x1 + 3*x2 + noise
    constexpr std::size_t NP = 2;
    ctrlpp::rls<double, NP> estimator;

    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.01);
    std::uniform_real_distribution<double> input(-1.0, 1.0);

    Eigen::Vector2d true_theta;
    true_theta << 2.0, 3.0;

    for (int i = 0; i < 200; ++i) {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = true_theta.dot(phi) + noise(gen);
        estimator.update(y, phi);
    }

    auto theta_hat = estimator.parameters();
    REQUIRE_THAT(theta_hat(0), WithinAbs(2.0, 0.1));
    REQUIRE_THAT(theta_hat(1), WithinAbs(3.0, 0.1));
}

TEST_CASE("RLS with forgetting tracks time-varying parameters") {
    constexpr std::size_t NP = 2;
    ctrlpp::rls_config<double, NP> cfg;
    cfg.lambda = 0.95;
    ctrlpp::rls<double, NP> estimator(cfg);

    std::mt19937 gen(123);
    std::normal_distribution<double> noise(0.0, 0.01);
    std::uniform_real_distribution<double> input(-1.0, 1.0);

    Eigen::Vector2d theta1;
    theta1 << 2.0, 3.0;
    Eigen::Vector2d theta2;
    theta2 << -1.0, 5.0;

    // Phase 1: true params = [2, 3]
    for (int i = 0; i < 100; ++i) {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = theta1.dot(phi) + noise(gen);
        estimator.update(y, phi);
    }

    // Phase 2: true params change to [-1, 5]
    for (int i = 0; i < 200; ++i) {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = theta2.dot(phi) + noise(gen);
        estimator.update(y, phi);
    }

    auto theta_hat = estimator.parameters();
    REQUIRE_THAT(theta_hat(0), WithinAbs(-1.0, 0.3));
    REQUIRE_THAT(theta_hat(1), WithinAbs(5.0, 0.3));
}

TEST_CASE("RLS covariance stays bounded under low excitation") {
    constexpr std::size_t NP = 2;
    constexpr double bound = 1e4;
    ctrlpp::rls_config<double, NP> cfg;
    cfg.cov_upper_bound = bound;
    ctrlpp::rls<double, NP> estimator(cfg);

    // Feed constant phi (zero excitation)
    Eigen::Vector2d phi;
    phi << 1.0, 0.0;

    for (int i = 0; i < 1000; ++i) {
        estimator.update(1.0, phi);
    }

    double trace = estimator.covariance().trace();
    REQUIRE(trace <= bound * static_cast<double>(NP) + 1e-6);
}

TEST_CASE("RLS covariance is symmetric after every update") {
    constexpr std::size_t NP = 3;
    ctrlpp::rls<double, NP> estimator;

    std::mt19937 gen(77);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (int i = 0; i < 50; ++i) {
        Eigen::Vector3d phi;
        phi << dist(gen), dist(gen), dist(gen);
        double y = dist(gen);
        estimator.update(y, phi);

        auto P = estimator.covariance();
        double asym = (P - P.transpose()).norm();
        REQUIRE(asym < 1e-12);
    }
}

TEST_CASE("RLS default construction produces reasonable initial state") {
    constexpr std::size_t NP = 2;
    ctrlpp::rls<double, NP> estimator;

    // Parameters should be zero-initialized
    REQUIRE_THAT(estimator.parameters()(0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(estimator.parameters()(1), WithinAbs(0.0, 1e-15));

    // Covariance should be large diagonal
    REQUIRE(estimator.covariance()(0, 0) > 100.0);
    REQUIRE(estimator.covariance()(1, 1) > 100.0);
}

TEST_CASE("RLS with forgetting factor 1.0 converges monotonically on stationary data") {
    constexpr std::size_t NP = 2;
    ctrlpp::rls_config<double, NP> cfg;
    cfg.lambda = 1.0;
    ctrlpp::rls<double, NP> estimator(cfg);

    std::mt19937 gen(99);
    std::normal_distribution<double> noise(0.0, 0.01);
    std::uniform_real_distribution<double> input(-1.0, 1.0);

    Eigen::Vector2d true_theta;
    true_theta << 1.5, -2.5;

    double prev_error = std::numeric_limits<double>::max();
    int monotonic_violations = 0;

    for (int i = 0; i < 100; ++i) {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = true_theta.dot(phi) + noise(gen);
        estimator.update(y, phi);

        double error = (estimator.parameters() - true_theta).norm();
        // Allow small noise-induced violations, but track them
        if (error > prev_error + 0.05) {
            ++monotonic_violations;
        }
        prev_error = error;
    }

    // With noise, a few violations are acceptable, but should be rare
    REQUIRE(monotonic_violations < 10);
}
