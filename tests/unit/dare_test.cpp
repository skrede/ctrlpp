#include "ctrlpp/dare.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>


TEST_CASE("dare scalar integrator golden ratio") {
    // A=1, B=1, Q=1, R=1 => P = (1+sqrt(5))/2 (golden ratio)
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    CHECK_THAT((*result)(0, 0), Catch::Matchers::WithinAbs(golden, 1e-10));
}

TEST_CASE("dare double integrator 2x2") {
    // Double integrator: A=[[1,1],[0,1]], B=[[0.5],[1]], Q=I, R=1
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0,
         0.0, 1.0;
    B << 0.5,
         1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto P = *result;

    // Verify P is symmetric
    CHECK((P - P.transpose()).norm() < 1e-10);

    // Verify P satisfies the DARE equation:
    // A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q = 0
    auto AtPA = A.transpose() * P * A;
    auto BtPB = B.transpose() * P * B;
    Eigen::Matrix<double, 1, 1> S = R + BtPB;
    auto AtPB = A.transpose() * P * B;
    auto residual = AtPA - P - AtPB * S.inverse() * AtPB.transpose() + Q;
    CHECK(residual.norm() < 1e-10);

    // Verify closed-loop stability: eigenvalues of A - B*K inside unit circle
    // K = (R + B^T P B)^{-1} B^T P A
    Eigen::Matrix<double, 1, 2> K = S.inverse() * B.transpose() * P * A;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for (int i = 0; i < 2; ++i)
        CHECK(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("dare 3-state system") {
    // 3-state system: A rotational + damping
    Eigen::Matrix<double, 3, 3> A, Q;
    Eigen::Matrix<double, 3, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.9, 0.1, 0.0,
         0.0, 0.8, 0.2,
         0.0, 0.0, 0.7;
    B << 0.0,
         0.0,
         1.0;
    Q = Eigen::Matrix<double, 3, 3>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 3, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto P = *result;

    // Verify symmetry
    CHECK((P - P.transpose()).norm() < 1e-10);

    // Verify DARE residual
    auto AtPA = A.transpose() * P * A;
    auto BtPB = B.transpose() * P * B;
    Eigen::Matrix<double, 1, 1> S = R + BtPB;
    auto AtPB = A.transpose() * P * B;
    auto residual = AtPA - P - AtPB * S.inverse() * AtPB.transpose() + Q;
    CHECK(residual.norm() < 1e-10);

    // Verify positive semi-definite (all eigenvalues >= 0)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> eigsolver(P);
    for (int i = 0; i < 3; ++i)
        CHECK(eigsolver.eigenvalues()(i) >= -1e-10);
}

TEST_CASE("dare non-stabilizable returns nullopt") {
    // A has unstable mode at eigenvalue 2, B cannot reach it
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 2.0, 0.0,
         0.0, 0.5;
    B << 0.0,
         1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK_FALSE(result.has_value());
}

TEST_CASE("dare with N cross-weight") {
    // Verify DARE with cross-weight N: solution should match transformed standard DARE
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B, N;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0,
         0.0, 1.0;
    B << 0.5,
         1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    N << 0.1,
         0.2;

    auto result_with_n = ctrlpp::dare<double, 2, 1>(A, B, Q, R, N);
    REQUIRE(result_with_n.has_value());

    // Transform manually: Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T
    auto Rinv = R.inverse();
    Eigen::Matrix<double, 2, 2> Qp = Q - N * Rinv * N.transpose();
    Eigen::Matrix<double, 2, 2> Ap = A - B * Rinv * N.transpose();

    auto result_standard = ctrlpp::dare<double, 2, 1>(Ap, B, Qp, R);
    REQUIRE(result_standard.has_value());

    CHECK((*result_with_n - *result_standard).norm() < 1e-10);
}

TEST_CASE("dare solution is symmetric") {
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.9, 0.1,
         0.0, 0.8;
    B << 1.0,
         0.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto P = *result;
    CHECK((P - P.transpose()).norm() < 1e-10);
}
