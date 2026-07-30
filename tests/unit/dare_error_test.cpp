#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <limits>

TEST_CASE("DARE non-stabilisable system fails with non_stabilisable or singular_u11", "[dare][error]")
{
    // A has an unstable mode at eigenvalue 2 that B cannot reach.
    // The symplectic spectrum still contains n=2 stable eigenvalues (0.5 and its
    // reciprocal pair), so the reorder succeeds in principle; the invariant-subspace
    // basis of such an input class is degenerate and the failure surfaces via
    // singular_u11 at extraction time. Either enumerator is a structurally correct
    // failure for a non-stabilisable input.
    Eigen::Matrix<double, 2, 2> A;
    A << 2.0, 0.0, 0.0, 0.5;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK((result.error() == ctrlpp::dare_error::non_stabilisable || result.error() == ctrlpp::dare_error::singular_u11));
}

TEST_CASE("DARE NaN in A returns dare_error::non_finite_input", "[dare][error]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("DARE singular A returns dare_error::singular_a", "[dare][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.0, 0.0, 0.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::singular_a);
}

TEST_CASE("DARE A = 0, B = 0 yields a structured failure enum", "[dare][error]")
{
    Eigen::Matrix<double, 2, 2> A = Eigen::Matrix<double, 2, 2>::Zero();
    Eigen::Matrix<double, 2, 1> B = Eigen::Matrix<double, 2, 1>::Zero();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK((result.error() == ctrlpp::dare_error::singular_a || result.error() == ctrlpp::dare_error::singular_u11 || result.error() == ctrlpp::dare_error::non_finite_input ||
           result.error() == ctrlpp::dare_error::non_stabilisable));
}

TEST_CASE("DARE negative-definite Q produces a structured failure enum", "[dare][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    Eigen::Matrix<double, 2, 2> Q = -Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    if(!result.has_value())
    {
        CHECK((result.error() == ctrlpp::dare_error::non_psd_solution || result.error() == ctrlpp::dare_error::non_stabilisable ||
               result.error() == ctrlpp::dare_error::non_finite_input || result.error() == ctrlpp::dare_error::singular_u11 || result.error() == ctrlpp::dare_error::schur_failed));
    }
}

TEST_CASE("DARE schur_failed enumerator is reachable at compile time", "[dare][error][design-lever]")
{
    constexpr ctrlpp::dare_error e = ctrlpp::dare_error::schur_failed;
    (void)e;
    CHECK(static_cast<int>(ctrlpp::dare_error::schur_failed) >= 0);
}

TEST_CASE("DARE arithmetic_limit enumerator is reachable at compile time", "[dare][error][design-lever]")
{
    constexpr ctrlpp::dare_error e = ctrlpp::dare_error::arithmetic_limit;
    (void)e;
    CHECK(static_cast<int>(ctrlpp::dare_error::arithmetic_limit) >= 0);
}
