#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"
#include "ctrlpp/control/lqr.h"

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

TEST_CASE("DARE refuses a singular R", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.0;

    auto const result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());

    // The symplectic build needs R^{-1} to form G = B R^{-1} B', exactly as it
    // needs A^{-T}, so a rank-deficient R gets the enumerator that names it
    // rather than one describing a symptom. Before that enumerator existed the
    // caller was told their input was non-finite, which is a false statement
    // about data they chose deliberately.
    CHECK(result.error() == ctrlpp::dare_error::singular_r);

    // The gain forwards it unchanged.
    auto const gain = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(gain.has_value());
    CHECK(gain.error() == ctrlpp::dare_error::singular_r);

    // The cross-weight overload inverts R before the symplectic build ever sees
    // it, so it carries its own copy of the test.
    Eigen::Matrix<double, 2, 1> N;
    N << 0.1, 0.2;
    auto const crossed = ctrlpp::dare<double, 2, 1>(A, B, Q, R, N);
    REQUIRE_FALSE(crossed.has_value());
    CHECK(crossed.error() == ctrlpp::dare_error::singular_r);
}

TEST_CASE("DARE refuses a rank-deficient R instead of solving a different problem",
          "[dare][hardening][negative]")
{
    // The dangerous half of the same defect, and the reason the test is a rank
    // test rather than a finiteness check. A rank-deficient but NONZERO R does
    // not make the QR solve produce infinities: it produces a least-squares
    // answer over the leading rank columns, which is finite. Without an explicit
    // rank test the build would form a G that is not B R^{-1} B', the solve
    // would run to completion, and the caller would be handed a confident
    // solution to a problem they did not pose.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> B;
    B << 0.5, 0.0, 1.0, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    auto R = Eigen::Matrix<double, 2, 2>::Zero().eval();
    R(0, 0) = 1.0;  // rank 1 of 2, and every entry finite

    REQUIRE(R.allFinite());

    auto const result = ctrlpp::dare<double, 2, 2>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::singular_r);
}

TEST_CASE("DARE accepts an R that is ill-conditioned but not singular",
          "[dare][hardening][robustness]")
{
    // The boundary the rank test must not overshoot. A weighting spanning ten
    // decades is a numerical-conditioning question, not a domain violation, and
    // refusing it would turn one into the other.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> B;
    B << 0.5, 0.0, 1.0, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    auto R = ctrlpp::test::ill_conditioned_2x2<double>(1e10);

    auto const result = ctrlpp::dare<double, 2, 2>(A, B, Q, R);
    REQUIRE(result.has_value());

    // Positive definiteness is asserted; the Riccati residual deliberately is
    // NOT, and the reason is recorded rather than the assertion quietly
    // omitted. Measured here: the residual is 1.34e-7 against a term scale of
    // 8.13, i.e. 1.6e-8 relative. That is seven orders above the
    // counted-operation budget the well-conditioned cases use, and it is not a
    // solver defect -- a forward error of about cond(R) * eps = 1e10 * 2.2e-16
    // = 2.2e-6 is what this conditioning buys, and the observed value sits two
    // decades INSIDE it. A counted-operation budget models rounding only and is
    // simply the wrong oracle in this regime; the right one is scaled by a
    // conditioning estimate, which has not been derived. Asserting the
    // counted-op budget here would fail on a correct solve, and widening it
    // until it passed would be fitting a constant to an observation.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(result->P);
    for(int i = 0; i < 2; ++i)
        CHECK(pes.eigenvalues()(i) > 0.0);
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

TEST_CASE("DARE solves an ill-conditioned but well-posed problem",
          "[dare][hardening][robustness]")
{
    // A = [[1,1],[0,1]] is controllable from B = [0.5; 1] (rank[B, AB] = 2), and
    // Q = diag(1, 1e-10) is positive definite -- barely -- which makes the pair
    // detectable. A unique stabilizing positive-definite solution therefore
    // exists, so a refusal here would be a well-posed problem reported
    // unsolvable, and finiteness alone would be far weaker than the property
    // the solution is supposed to have.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto const result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto const& P = result->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(P);
    for(int i = 0; i < 2; ++i)
        CHECK(pes.eigenvalues()(i) > 0.0);

    auto const K = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, P);
    Eigen::Matrix<double, 2, 2> Acl = (A - B * K).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> ces(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(ces.eigenvalues()(i)) < 1.0);
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
