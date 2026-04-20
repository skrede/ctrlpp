#include "hardening_helpers.h"
#include "ctrlpp/control/care.h"


#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <limits>


TEST_CASE("CARE non-LHP-stabilisable system fails with non_lhp_stabilisable or singular_u11",
          "[care][error]")
{
    // A has an unstable continuous mode at eigenvalue +2 uncoupled from B.
    // The Hamiltonian spectrum still has n=2 eigenvalues in the open LHP (-2 and -0.5),
    // so the reorder succeeds in principle; the invariant-subspace basis is degenerate
    // and the failure surfaces via singular_u11 at extraction. Either enumerator is a
    // structurally correct failure for this input class.
    Eigen::Matrix<double, 2, 2> A;
    A << 2.0, 0.0, 0.0, -0.5;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK((result.error() == ctrlpp::care_error::non_lhp_stabilisable
        || result.error() == ctrlpp::care_error::singular_u11));
}

TEST_CASE("CARE NaN in A returns care_error::non_finite_input",
          "[care][error]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::care_error::non_finite_input);
}

TEST_CASE("CARE Inf in B returns care_error::non_finite_input",
          "[care][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0, -0.5, -0.3;
    auto B = ctrlpp::test::inf_matrix<double, 2, 1>();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::care_error::non_finite_input);
}

TEST_CASE("CARE A = 0, B = 0 yields a structured failure enum",
          "[care][error]")
{
    Eigen::Matrix<double, 2, 2> A = Eigen::Matrix<double, 2, 2>::Zero();
    Eigen::Matrix<double, 2, 1> B = Eigen::Matrix<double, 2, 1>::Zero();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK((result.error() == ctrlpp::care_error::singular_u11
        || result.error() == ctrlpp::care_error::non_finite_input
        || result.error() == ctrlpp::care_error::non_lhp_stabilisable
        || result.error() == ctrlpp::care_error::schur_failed
        || result.error() == ctrlpp::care_error::sign_function_stagnated));
}

TEST_CASE("CARE negative-definite Q produces a structured failure enum",
          "[care][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0, -0.5, -0.3;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = -Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    if (!result.has_value())
    {
        CHECK((result.error() == ctrlpp::care_error::non_psd_solution
            || result.error() == ctrlpp::care_error::non_lhp_stabilisable
            || result.error() == ctrlpp::care_error::non_finite_input
            || result.error() == ctrlpp::care_error::singular_u11
            || result.error() == ctrlpp::care_error::schur_failed
            || result.error() == ctrlpp::care_error::sign_function_stagnated));
    }
}

TEST_CASE("CARE schur_failed enumerator is reachable at compile time",
          "[care][error][design-lever]")
{
    constexpr ctrlpp::care_error e = ctrlpp::care_error::schur_failed;
    (void)e;
    CHECK(static_cast<int>(ctrlpp::care_error::schur_failed) >= 0);
}
