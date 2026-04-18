// Correctness and no-allocation hardening for ctrlpp::lqr_gain_continuous.
//
// lqr_gain_continuous(A, B, Q, R) computes K = R^{-1} B^T P where P solves the
// continuous-time Riccati equation A^T P + P A - P B R^{-1} B^T P + Q = 0. The
// closed-loop dynamics are A - B K, and for a stabilisable pair (A, B) all closed-
// loop eigenvalues must lie in the open left half-plane.

#define EIGEN_RUNTIME_NO_MALLOC

#include "ctrlpp/control/lqr.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


TEST_CASE("lqr_gain_continuous stabilises the continuous double integrator",
          "[lqr][continuous]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0, 0.0, 0.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto K = ctrlpp::lqr_gain_continuous<double, 2, 1>(A, B, Q, R);
    REQUIRE(K.has_value());

    Eigen::Matrix<double, 2, 2> Acl = A - B * (*K);
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> eig(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(eig.eigenvalues()(i).real() < 0.0);
}

TEST_CASE("lqr_gain_continuous matches analytic gain on scalar system",
          "[lqr][continuous]")
{
    // A=0, B=1, Q=1, R=1 -> P^2 = 1 -> P = 1, K = R^{-1} B^T P = 1
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 0.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto K = ctrlpp::lqr_gain_continuous<double, 1, 1>(A, B, Q, R);
    REQUIRE(K.has_value());
    CHECK_THAT((*K)(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("lqr_gain_continuous returns nullopt on non-LHP-stabilisable system",
          "[lqr][continuous][negative]")
{
    // A has an unstable mode at +2 uncoupled from B.
    Eigen::Matrix<double, 2, 2> A;
    A << 2.0, 0.0, 0.0, -0.5;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto K = ctrlpp::lqr_gain_continuous<double, 2, 1>(A, B, Q, R);
    CHECK(!K.has_value());
}

TEST_CASE("lqr_gain_continuous hot path performs zero heap allocation (NX=4, NU=2)",
          "[lqr][continuous][hardening][nomalloc]")
{
    Eigen::Matrix<double, 4, 4> A = Eigen::Matrix<double, 4, 4>::Zero();
    for(int i = 0; i < 4; ++i)
        A(i, i) = -0.5;
    for(int i = 0; i + 1 < 4; ++i)
        A(i, i + 1) = 1.0;

    Eigen::Matrix<double, 4, 2> B = Eigen::Matrix<double, 4, 2>::Zero();
    B(1, 0) = 1.0;
    B(3, 1) = 1.0;

    Eigen::Matrix<double, 4, 4> Q = Eigen::Matrix<double, 4, 4>::Identity();
    Eigen::Matrix<double, 2, 2> R = 0.1 * Eigen::Matrix<double, 2, 2>::Identity();

    auto warmup = ctrlpp::lqr_gain_continuous<double, 4, 2>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto K = ctrlpp::lqr_gain_continuous<double, 4, 2>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(K.has_value());
}
