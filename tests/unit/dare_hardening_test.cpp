// What the oracles in this file decide.
//
// The solver returns a matrix that is supposed to SOLVE an equation, so that is
// what is asserted:
//
//  * Every accepted solution is held to the Riccati residual against a
//    counted-operation budget scaled by the largest of the four terms that
//    cancel to produce it. Positive definiteness alone does not identify the
//    solution -- a positive definite matrix that solves nothing passes it -- so
//    definiteness is asserted alongside the residual, never instead of it.
//  * Positive definiteness is asserted with the exact contract boundary. The
//    floor is a floor at zero, not below it: slack in the direction that admits
//    a negative eigenvalue admits the very matrix the case is named against.
//  * Symmetry is asserted EXACTLY. The solver symmetrizes the raw quotient
//    before returning it, so a bitwise symmetric matrix is the contract and a
//    tolerance would admit one that is not.
//  * Every refusal names its enumerator. "No value" alone does not say the
//    solver diagnosed the caller's actual fault.
//
// What they deliberately do not decide. One case here accepts an ill-conditioned
// weighting and does NOT assert the residual; the reason is recorded at that
// case rather than the assertion quietly omitted, and it is that a
// counted-operation budget models rounding only and is the wrong oracle once
// the conditioning dominates.

#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"
#include "ctrlpp/control/lqr.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

TEST_CASE("DARE refuses a non-finite state matrix", "[dare][hardening][negative]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    // The old name promised "nullopt or NaN" while the assertion below demanded
    // the first alternative unconditionally, so the name described a weaker
    // contract than the test enforced. The enumerator is asserted too: a refusal
    // that does not say WHY sends the caller looking in the wrong place.
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("DARE refuses a non-finite input matrix", "[dare][hardening][negative]")
{
    auto A = Eigen::Matrix<double, 2, 2>::Identity();
    auto B = ctrlpp::test::nan_matrix<double, 2, 1>();
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
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

    auto const& P = result->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    // The case proved P was positive definite and never that P solves anything,
    // so any positive definite matrix of the right size passed it. The residual
    // is the property that identifies the solution.
    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    // Verify positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(P);
    for(int i = 0; i < 2; ++i)
        CHECK(eigsolver.eigenvalues()(i) > 0.0);

    // Symmetry is exact, not a tolerance: the solver symmetrizes the raw
    // quotient U21 * U11^-1 before returning it, so the two triangles hold the
    // same bits. An unexplained 1e-14 admitted an asymmetry the construction
    // cannot produce and would have hidden a dropped symmetrization.
    REQUIRE((P - P.transpose()).norm() == 0.0);
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

    constexpr double eps = std::numeric_limits<double>::epsilon();

    // Analytical: P = golden ratio = (1+sqrt(5))/2. The budget is carried from
    // the residual to the solution through the residual's own derivative: for
    // this data r(P) = 1 - P^2 / (1 + P), so r'(P) = -(P^2 + 2P) / (1 + P)^2,
    // which is 0.854 at the golden ratio. A residual inside the counted chain
    // therefore puts the solution inside that chain divided by 0.854. Two
    // further roundings enter on the test side, the square root and the sum.
    double const golden = (1.0 + std::sqrt(5.0)) / 2.0;
    double const residual_slope = (golden * golden + 2.0 * golden)
                                  / ((1.0 + golden) * (1.0 + golden));
    constexpr int analytic_ops = 2;
    double const budget =
        ctrlpp::test::riccati_residual_ops<1, 1> * eps * golden / residual_slope
        + analytic_ops * eps * golden;

    CAPTURE(result->P(0, 0), golden, budget);
    REQUIRE(std::abs(result->P(0, 0) - golden) <= budget);
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

    auto const& P = result->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    // The floor is at zero, where the contract is. The previous form admitted an
    // eigenvalue down to -1e-10 -- half a million units in the last place of
    // slack pointing INTO the indefinite half-space -- so a solution that was
    // not positive definite passed a case named for positive definiteness. Both
    // eigenvalues here are of order one, nowhere near the boundary, so nothing
    // is being tightened onto a knife edge.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(P);
    for(int i = 0; i < 2; ++i)
    {
        CAPTURE(i, eigsolver.eigenvalues()(i));
        CHECK(eigsolver.eigenvalues()(i) > 0.0);
    }

    REQUIRE((P - P.transpose()).norm() == 0.0);
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

TEST_CASE("DARE refuses a non-finite state weighting", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}
