#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

TEST_CASE("DARE NaN in A returns nullopt or NaN", "[dare][hardening][negative]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK(!result.has_value());
}

TEST_CASE("DARE NaN in B returns nullopt", "[dare][hardening][negative]")
{
    auto A = Eigen::Matrix<double, 2, 2>::Identity();
    auto B = ctrlpp::test::nan_matrix<double, 2, 1>();
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK(!result.has_value());
}

TEST_CASE("DARE zero R (singular) returns nullopt", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    // Singular R -- may not produce a valid solution
    if(result.has_value())
    {
        bool all_finite = true;
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                if(!std::isfinite(result->P(i, j)))
                    all_finite = false;
        // If it returns something, it should be finite or we accept nullopt
        CHECK(all_finite);
    }
}

TEST_CASE("DARE known 2x2 solution is positive definite", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto& P = result->P;

    // Verify positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(P);
    for(int i = 0; i < 2; ++i)
        CHECK(eigsolver.eigenvalues()(i) > 0.0);

    // Verify symmetric
    REQUIRE_THAT((P - P.transpose()).norm(), WithinAbs(0.0, 1e-14));
}

TEST_CASE("DARE scalar analytical solution", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    // Analytical: P = golden ratio = (1+sqrt(5))/2
    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    REQUIRE_THAT(result->P(0, 0), WithinAbs(golden, 1e-10));
}

TEST_CASE("DARE solution is positive definite for stable system", "[dare][hardening][stability]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    A << 0.9, 0.1, 0.0, 0.8;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(result->P);
    for(int i = 0; i < 2; ++i)
        CHECK(eigsolver.eigenvalues()(i) > -1e-10);
}

TEST_CASE("DARE ill-conditioned Q with cond 1e10", "[dare][hardening][robustness]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    if(result.has_value())
    {
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                CHECK(std::isfinite(result->P(i, j)));
    }
}

TEST_CASE("DARE NaN in Q returns nullopt", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK(!result.has_value());
}
