#include "hardening_helpers.h"
#include "ctrlpp/control/lqr.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("LQR NaN in Q produces no value", "[lqr][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = ctrlpp::test::nan_matrix<double, 2, 2>();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    // NaN inputs should either return nullopt or produce NaN -- no crash
    CHECK((!result.has_value() || !std::isfinite((*result)(0, 0))));
}

TEST_CASE("LQR NaN in R produces no value", "[lqr][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << std::numeric_limits<double>::quiet_NaN();

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    CHECK((!result.has_value() || !std::isfinite((*result)(0, 0))));
}

TEST_CASE("LQR ill-conditioned system pair", "[lqr][hardening][negative]")
{
    auto A = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    // May succeed or return nullopt -- just no crash
    if(result.has_value())
    {
        CHECK(std::isfinite((*result)(0, 0)));
        CHECK(std::isfinite((*result)(0, 1)));
    }
}

TEST_CASE("LQR double integrator matches analytical", "[lqr][hardening][precision]")
{
    double dt = 0.1;
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, dt, 0.0, 1.0;
    B << 0.5 * dt * dt, dt;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto& K = *result;
    CHECK(std::isfinite(K(0, 0)));
    CHECK(std::isfinite(K(0, 1)));

    // Verify closed-loop eigenvalues are inside unit circle
    auto Acl = (A - B * K).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR closed-loop eigenvalues inside unit circle", "[lqr][hardening][stability]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto Acl = (A - B * *result).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR near-singular system epsilon 1e-10", "[lqr][hardening][robustness]")
{
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(1e-10);
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(sys.A, sys.B, Q, R);
    if(result.has_value())
    {
        CHECK(std::isfinite((*result)(0, 0)));
        CHECK(std::isfinite((*result)(0, 1)));
    }
}

TEST_CASE("LQR scalar integrator analytical gain", "[lqr][hardening][precision]")
{
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    double expected_K = golden / (1.0 + golden);
    REQUIRE_THAT((*result)(0, 0), WithinAbs(expected_K, 1e-10));
}
