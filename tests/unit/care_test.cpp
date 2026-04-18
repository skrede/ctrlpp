#include "ctrlpp/control/care.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>


TEST_CASE("care scalar integrator analytical")
{
    // A=0, B=1, Q=1, R=1 => P^2 = 1 => P = 1
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 0.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    CHECK_THAT((*result)(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("care double integrator analytical")
{
    // Continuous double integrator: A=[[0,1],[0,0]], B=[[0],[1]], Q=I, R=1
    // Analytical P = [[sqrt(3), 1], [1, sqrt(3)]]
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.0, 1.0, 0.0, 0.0;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    const auto& P = *result;
    const double s3 = std::sqrt(3.0);

    CHECK_THAT(P(0, 0), Catch::Matchers::WithinAbs(s3, 1e-10));
    CHECK_THAT(P(1, 1), Catch::Matchers::WithinAbs(s3, 1e-10));
    CHECK_THAT(P(0, 1), Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(P(1, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));

    auto residual = A.transpose() * P + P * A - P * B * R.inverse() * B.transpose() * P + Q;
    CHECK(residual.norm() < 1e-10);
}

TEST_CASE("care 3-state damped system")
{
    Eigen::Matrix<double, 3, 3> A, Q;
    Eigen::Matrix<double, 3, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << -0.5, 1.0, 0.0,
         0.0, -0.8, 1.0,
         0.0, 0.0, -1.0;
    B << 0.0, 0.0, 1.0;
    Q = Eigen::Matrix<double, 3, 3>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 3, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    const auto& P = *result;

    CHECK((P - P.transpose()).norm() < 1e-10);

    auto residual = A.transpose() * P + P * A - P * B * R.inverse() * B.transpose() * P + Q;
    CHECK(residual.norm() < 1e-10);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> eigsolver(P);
    for(int i = 0; i < 3; ++i)
        CHECK(eigsolver.eigenvalues()(i) >= -1e-10);
}

TEST_CASE("care closed-loop stable")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.0, 1.0, -0.1, -0.2;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    const auto& P = *result;
    Eigen::Matrix<double, 1, 2> K = R.inverse() * B.transpose() * P;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;

    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(solver.eigenvalues()(i).real() < 0.0);
}

TEST_CASE("care non-stabilizable returns nullopt")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 2.0, 0.0, 0.0, -0.5;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    CHECK_FALSE(result.has_value());
}

TEST_CASE("care with N cross-weight")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B, N;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.0, 1.0, -0.5, -0.3;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    N << 0.1, 0.2;

    auto result_with_n = ctrlpp::care<double, 2, 1>(A, B, Q, R, N);
    REQUIRE(result_with_n.has_value());

    auto Rinv = R.inverse();
    Eigen::Matrix<double, 2, 2> Qp = Q - N * Rinv * N.transpose();
    Eigen::Matrix<double, 2, 2> Ap = A - B * Rinv * N.transpose();

    auto result_standard = ctrlpp::care<double, 2, 1>(Ap, B, Qp, R);
    REQUIRE(result_standard.has_value());

    CHECK((*result_with_n - *result_standard).norm() < 1e-10);
}

TEST_CASE("care solution is symmetric")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << -0.5, 1.0, 0.0, -0.8;
    B << 1.0, 0.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    CHECK((*result - result->transpose()).norm() < 1e-10);
}
