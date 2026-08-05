// What the oracles in this file decide.
//
// The solver returns a matrix that is supposed to SOLVE an equation, so that is
// what is asserted:
//
//  * The ACCEPTANCE contract is a statement about the ANSWER: an accepted
//    solution retains more than half of the scalar type's significand. It is
//    asserted through the library's own shared rule and never re-spelled, so an
//    anchor cannot pass against its own copy of the contract while the contract
//    itself has moved.
//  * Well-conditioned anchors are ALSO held to the Riccati residual against a
//    counted-operation budget scaled by the largest of the four terms that
//    cancel to produce it. That is a tighter, rounding-floor statement that
//    those particular poses satisfy; it is not the acceptance contract, and it
//    is not asserted on poses whose conditioning dominates their rounding.
//    Positive definiteness alone does not identify the solution -- a positive
//    definite matrix that solves nothing passes it -- so definiteness is
//    asserted alongside, never instead.
//  * Positive definiteness is asserted with the exact contract boundary. The
//    floor is a floor at zero, not below it: slack in the direction that admits
//    a negative eigenvalue admits the very matrix the case is named against.
//  * Symmetry is asserted EXACTLY. The solver symmetrizes the raw quotient
//    before returning it, so a bitwise symmetric matrix is the contract and a
//    tolerance would admit one that is not.
//  * Every refusal names its enumerator. "No value" alone does not say the
//    solver diagnosed the caller's actual fault.
//
// What they deliberately do not decide. No case here asserts a residual bound as
// if it were an accuracy bound. On the population this solver is measured
// against, the answers that keep half the significand and the answers that do
// not are contiguous on the residual, so a threshold on that quantity separates
// nothing; the case that reaches an ill-conditioned weighting asserts the
// accuracy rule and the refusal boundary instead, with the measured forward
// errors recorded at the case.

#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"
#include "ctrlpp/control/lqr.h"
#include "ctrlpp/detail/riccati_solution.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>

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

TEST_CASE("DARE cross-weight overload classifies every non-finite operand", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Matrix<double, 2, 1> N;
    N << 0.1, 0.2;

    auto require_non_finite = [](auto const &a, auto const &b, auto const &q, auto const &r, auto const &n)
    {
        auto const result = ctrlpp::dare<double, 2, 1>(a, b, q, r, n);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
    };

    auto bad_A  = A;
    bad_A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require_non_finite(bad_A, B, Q, R, N);

    auto bad_B  = B;
    bad_B(0, 0) = std::numeric_limits<double>::infinity();
    require_non_finite(A, bad_B, Q, R, N);

    auto bad_Q  = Q;
    bad_Q(0, 0) = -std::numeric_limits<double>::infinity();
    require_non_finite(A, B, bad_Q, R, N);

    auto bad_R  = R;
    bad_R(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require_non_finite(A, B, Q, bad_R, N);

    auto bad_N  = N;
    bad_N(0, 0) = std::numeric_limits<double>::infinity();
    require_non_finite(A, B, Q, R, bad_N);
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

TEST_CASE("DARE refuses a rank-deficient R instead of solving a different problem", "[dare][hardening][negative]")
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
    auto Q  = Eigen::Matrix<double, 2, 2>::Identity();
    auto R  = Eigen::Matrix<double, 2, 2>::Zero().eval();
    R(0, 0) = 1.0; // rank 1 of 2, and every entry finite

    REQUIRE(R.allFinite());

    auto const result = ctrlpp::dare<double, 2, 2>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::singular_r);
}

TEST_CASE("DARE separates an ill-conditioned R from a rank-deficient one", "[dare][hardening][robustness]")
{
    // The boundary the rank test must not overshoot. A weighting spanning many
    // decades is a numerical-conditioning question, not a domain violation, and
    // reporting it as a rank deficiency would turn one into the other. That is
    // what this case guards, and it is asserted on EVERY conditioning below --
    // no outcome here is ever `singular_r`.
    //
    // What the outcome is instead is a precision question with a measured
    // answer. The accepted set on this family ends where binary64 stops being
    // able to deliver half a significand, and both sides of that edge are
    // asserted rather than one. Relative forward errors against an independent
    // extended-precision solution of the identical pose, computed by
    // Newton-Kleinman rather than by a Schur decomposition:
    //
    //     cond(R)   forward error     half-significand line is 1.4901e-08
    //     1e2       1.1778e-14
    //     1e6       1.1855e-10
    //     1e7       4.5861e-10        last conditioning that keeps half
    //     1e8       1.5757e-08        first that does not
    //     1e9       6.4807e-08
    //     1e10      4.5108e-08
    //
    // The solver used to return the last three. They are not correct answers
    // that were lost: every one of them has lost more than half of binary64's
    // significand, by an independent reference, and the estimate the solver now
    // forms agrees with that reference to five significant figures on all six
    // rows. Above 1e10 the pose was already refused before this edge existed.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> B;
    B << 0.5, 0.0, 1.0, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();

    int accepted = 0;
    int declined = 0;
    for(const double conditioning : {1e2, 1e4, 1e6, 1e7, 1e8, 1e9, 1e10, 1e12})
    {
        CAPTURE(conditioning);
        auto R            = ctrlpp::test::ill_conditioned_2x2<double>(conditioning);
        auto const result = ctrlpp::dare<double, 2, 2>(A, B, Q, R);

        if(result.has_value())
        {
            ++accepted;
            CHECK(conditioning <= 1e7);
            // An accepted answer keeps more than half the significand, by the
            // one shared rule, and is positive definite.
            CHECK(ctrlpp::test::riccati_accuracy_of<double, 2, 2>(A, B, Q, R, result->P)
                  == ctrlpp::detail::riccati_accuracy::within);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(result->P);
            for(int i = 0; i < 2; ++i)
                CHECK(pes.eigenvalues()(i) > 0.0);
        }
        else
        {
            ++declined;
            CHECK(conditioning >= 1e8);
            // The point of the case: never a rank deficiency.
            CHECK(result.error() != ctrlpp::dare_error::singular_r);
            CHECK(result.error() == ctrlpp::dare_error::arithmetic_limit);
        }
    }
    CHECK(accepted > 0);
    CHECK(declined > 0);
}

TEST_CASE("DARE names the weight ratio that dissolves the input weighting", "[dare][hardening][precision]")
{
    // The solve factorizes the symplectic operands ONCE, from the equilibrated
    // weighting, because that is the operand the symplectic build consumes. A
    // consequence reaches the caller: the rank verdict on the input weighting is
    // made at the equilibrated scale, so a weighting that is nonzero on its own
    // but vanishes against the state weighting is reported as what it is.
    //
    // These poses used to report `arithmetic_limit`, which sends a caller to look
    // at precision when the obstacle is the weight ratio they chose. Both halves
    // are asserted: the ratio that dissolves the weighting names it, and a ratio
    // wide enough to be uncomfortable but not wide enough to dissolve it is still
    // solved -- without the second the case would pass on an implementation that
    // called every wide ratio singular.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.1, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.005, 0.1;

    auto weighting = [](double q_scale, double r_scale)
    {
        Eigen::Matrix<double, 2, 2> Q = q_scale * Eigen::Matrix<double, 2, 2>::Identity();
        Eigen::Matrix<double, 1, 1> R;
        R << r_scale;
        return std::pair{Q, R};
    };

    std::size_t dissolved{};
    for(const double q_scale : {1e150, 1e250, 1e300, 1e308})
    {
        for(const double r_scale : {1e-250, 1e-300, 1e-320})
        {
            CAPTURE(q_scale, r_scale);
            // The ratio exceeds the type's range, so the equilibrated weighting
            // underflows to zero whatever the pose.
            REQUIRE(r_scale / q_scale == 0.0);
            const auto [Q, R] = weighting(q_scale, r_scale);
            const auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
            REQUIRE_FALSE(result.has_value());
            CHECK(result.error() == ctrlpp::dare_error::singular_r);
            ++dissolved;
        }
    }
    CHECK(dissolved == 12);

    // A ratio of eight decades: wide, and solved. Nothing here is refused for
    // being merely lopsided, and the pose family's own accepted ceiling sits at
    // this ratio -- measured, and unchanged by the reordering above.
    const auto [Q_wide, R_wide] = weighting(1e4, 1e-4);
    const auto wide             = ctrlpp::dare<double, 2, 1>(A, B, Q_wide, R_wide);
    REQUIRE(wide.has_value());
    CHECK(wide->K.allFinite());
}

TEST_CASE("DARE returns the gain it verified rather than one re-formed from P", "[dare][hardening][precision]")
{
    // WHAT THIS PINS IS THE DIFFERENCE BETWEEN THE TWO, ON THE POSE WHERE IT IS
    // TOTAL RATHER THAN IN THE LAST BITS.
    //
    // The gain is homogeneous of degree zero in (P, Q, R), so the equilibrated
    // gain the solve verified IS the caller's gain. Re-forming it from the
    // returned P at the caller's own scale is a different computation with a
    // different failure mode: at the top of the range `R + B'PB` leaves the range
    // while P itself is an ordinary normal number, and a rank-revealing solve
    // handed that sum returns a ZERO gain. Zero gain is not a refusal and not a
    // small error -- it is no feedback at all, returned as if it were the control
    // law.
    //
    // The scalar pose makes the arithmetic checkable by hand: for A = 0.5, B = 1
    // and equal weights the scale-free gain is 0.2655644370746374, and the
    // returned P at common scale c is c times the unit-scale P.
    Eigen::Matrix<double, 1, 1> A;
    A << 0.5;
    Eigen::Matrix<double, 1, 1> B;
    B << 1.0;
    Eigen::Matrix<double, 1, 1> weight;
    weight << 1e308;

    const auto solved = ctrlpp::dare<double, 1, 1>(A, B, weight, weight);
    REQUIRE(solved.has_value());

    // The re-forming a caller would write, and what it produces here.
    const auto BtP     = (B.transpose() * solved->P).eval();
    const auto sum     = (weight + BtP * B).eval();
    const auto reformed = sum.colPivHouseholderQr().solve(BtP * A).eval();
    CHECK(sum(0, 0) == std::numeric_limits<double>::infinity());
    CHECK(reformed(0, 0) == 0.0);

    // The solve's own gain on the same pose, against the analytical value.
    const double analytic = 0.2655644370746374;
    CHECK(solved->K.allFinite());
    CHECK(std::abs(solved->K(0, 0) - analytic) <= 4.0 * std::numeric_limits<double>::epsilon() * analytic);

    // And the public gain helper hands back the solve's, not the re-forming.
    const auto helper = ctrlpp::lqr_gain<double, 1, 1>(A, B, weight, weight);
    REQUIRE(helper.has_value());
    CHECK((*helper)(0, 0) == solved->K(0, 0));
}

TEST_CASE("DARE accuracy helper is not blind to exchanged weightings", "[dare][hardening][precision]")
{
    // WHAT IS PINNED HERE IS AN ARGUMENT ORDER, BY BEHAVIOR RATHER THAN BY
    // CONVENTION.
    //
    // The accuracy helper used to take its two weightings input-first while its
    // own file's residual helper, the library's verification routine and the
    // public solver all take them state-first. At a square instantiation -- which
    // is the instantiation that exists, two states and two inputs -- the two
    // weighting types are the same type, so exchanging them at a call site
    // compiles, runs, and reports the accuracy of a DIFFERENT problem than the
    // one the case posed. Reordering the declaration alone would be a rename: it
    // makes the hazard less likely to be written, not detectable once written.
    //
    // So the order is asserted. The helper is handed a pose the solver accepted,
    // and it must say the answer keeps more than half the significand; handed the
    // same pose with the two weightings exchanged, it must NOT, because the
    // solution to one problem is not a solution to the other. Both directions are
    // required: without the first the case could pass on a helper that never
    // returns `within`, and without the second it could pass on a helper that
    // ignores its weightings entirely.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.1, 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> B;
    B << 0.005, 0.0, 0.1, 1.0;

    // Separated by four decades so the exchange is not a perturbation of the
    // posed problem but a different one, and so the case does not turn on where
    // a margin happens to fall.
    const auto Q = (1e4 * Eigen::Matrix<double, 2, 2>::Identity()).eval();
    const auto R = Eigen::Matrix<double, 2, 2>::Identity().eval();

    const auto solved = ctrlpp::dare<double, 2, 2>(A, B, Q, R);
    REQUIRE(solved.has_value());

    CHECK(ctrlpp::test::riccati_accuracy_of<double, 2, 2>(A, B, Q, R, solved->P)
          == ctrlpp::detail::riccati_accuracy::within);
    CHECK(ctrlpp::test::riccati_accuracy_of<double, 2, 2>(A, B, R, Q, solved->P)
          != ctrlpp::detail::riccati_accuracy::within);
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

    auto const &P        = result->P;
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
    double const golden         = (1.0 + std::sqrt(5.0)) / 2.0;
    double const residual_slope = (golden * golden + 2.0 * golden) / ((1.0 + golden) * (1.0 + golden));
    constexpr int analytic_ops  = 2;
    double const budget         = ctrlpp::test::riccati_residual_ops<1, 1> * eps * golden / residual_slope + analytic_ops * eps * golden;

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

    auto const &P        = result->P;
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

TEST_CASE("DARE solves an ill-conditioned but well-posed problem", "[dare][hardening][robustness]")
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

    auto const &P        = result->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(P);
    for(int i = 0; i < 2; ++i)
        CHECK(pes.eigenvalues()(i) > 0.0);

    auto const K                    = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, P);
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

TEST_CASE("DARE gain is invariant under common weight scaling", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.1, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.005, 0.1;

    const double gain_relative_margin = ctrlpp::test::riccati_gain_agreement_margin<double, 2, 1>();
    std::size_t accepted{};
    std::size_t arithmetic_declines{};
    std::size_t band_samples{};

    auto compare_common_scale = [&](const Eigen::Matrix<double, 2, 2> &Q, const Eigen::Matrix<double, 1, 1> &R, double common_scale)
    {
        const auto direct = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
        const auto scaled = ctrlpp::dare<double, 2, 1>(A, B, (Q / common_scale).eval(), (R / common_scale).eval());

        if(!direct || !scaled)
        {
            REQUIRE_FALSE(direct.has_value());
            REQUIRE_FALSE(scaled.has_value());
            CHECK(direct.error() == scaled.error());
            // This equality is a PREDICATE, not a check: it classifies a
            // decline rather than asserting one, so it decides which branch
            // runs. It is kept as an equality deliberately, and it is not
            // widened, because the census it feeds is asserted only to be
            // non-empty. Which individual points reach the accuracy enumerator
            // does move with the arithmetic -- six of the forty-two swept
            // points on some configurations and seven on others, measured
            // across four compilers, four optimization levels, fused
            // multiply-add off and on, and two releases of the linear-algebra
            // library -- but the count is at least six in every one of those
            // sixty-four, so the assertion below it holds on all of them. The
            // two posings agreeing on the enumerator, checked above, holds in
            // all sixty-four as well.
            if(direct.error() == ctrlpp::dare_error::arithmetic_limit)
                ++arithmetic_declines;
            return;
        }

        ++accepted;
        const auto direct_gain  = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, direct->P);
        const auto scaled_gain  = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, (R / common_scale).eval(), scaled->P);
        const double gain_scale = std::max(direct_gain.norm(), scaled_gain.norm());
        CAPTURE(common_scale, direct_gain, scaled_gain, gain_scale);
        CHECK((direct_gain - scaled_gain).norm() <= gain_relative_margin * gain_scale);
    };

    for(const double exponent : {-300.0, -200.0, -100.0, -18.0, -12.0, -6.0, 0.0, 2.0, 4.0, 6.0, 8.0, 8.2, 8.4, 8.6, 9.0, 10.0, 12.0, 18.0, 50.0, 100.0, 200.0, 300.0})
    {
        const double state_scale      = std::pow(10.0, exponent);
        Eigen::Matrix<double, 2, 2> Q = state_scale * Eigen::Matrix<double, 2, 2>::Identity();
        Eigen::Matrix<double, 1, 1> R;
        R << 0.1;
        compare_common_scale(Q, R, state_scale);
        if(exponent >= 4.0)
            ++band_samples;
    }

    for(const double exponent : {-300.0, -200.0, -100.0, -18.0, -12.0, -6.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 50.0, 100.0, 200.0, 300.0})
    {
        const double input_scale      = 0.1 * std::pow(10.0, exponent);
        Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
        Eigen::Matrix<double, 1, 1> R;
        R << input_scale;
        compare_common_scale(Q, R, input_scale);
        if(exponent >= 6.0)
            ++band_samples;
    }

    CHECK(band_samples > 0);
    CHECK(accepted > 0);
    CHECK(arithmetic_declines > 0);
}

TEST_CASE("DARE repairs or declines the measured state-heavy cases", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.1, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.005, 0.1;
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;

    const auto solved    = ctrlpp::dare<double, 2, 1>(A, B, (1e8 * Eigen::Matrix<double, 2, 2>::Identity()).eval(), R);
    const auto reference = ctrlpp::dare<double, 2, 1>(A, B, Eigen::Matrix<double, 2, 2>::Identity(), (R / 1e8).eval());
    REQUIRE(solved.has_value());
    REQUIRE(reference.has_value());

    const auto solved_gain            = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, solved->P);
    const auto reference_gain         = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, (R / 1e8).eval(), reference->P);
    const double gain_relative_margin = ctrlpp::test::riccati_gain_agreement_margin<double, 2, 1>();
    CHECK((solved_gain - reference_gain).norm() <= gain_relative_margin * reference_gain.norm());

    const auto beyond_precision = ctrlpp::dare<double, 2, 1>(A, B, (1e10 * Eigen::Matrix<double, 2, 2>::Identity()).eval(), R);
    REQUIRE_FALSE(beyond_precision.has_value());
    CHECK(beyond_precision.error() == ctrlpp::dare_error::arithmetic_limit);
}

TEST_CASE("DARE preserves comfortable common-scaled problems", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 0.8, 0.1, 0.0, 0.7;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.2, 0.4;

    std::size_t accepted{};
    for(int exponent = -12; exponent <= 7; ++exponent)
    {
        const double scale = std::pow(10.0, static_cast<double>(exponent));
        const auto Q       = (scale * Eigen::Matrix<double, 2, 2>::Identity()).eval();
        Eigen::Matrix<double, 1, 1> R;
        R << scale;
        const auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
        CAPTURE(exponent, scale);
        REQUIRE(result.has_value());
        ++accepted;
    }
    CHECK(accepted == 20);
}

namespace
{

// The scalar pose the common-scale claim is argued on: A = 0.5, B = 1, both
// weightings at a common positive scale. The exact solution is linear in that
// scale, so the unit-scale solve is an exact oracle for every other scale, and
// the gain it implies is the same at every scale -- which is the identity
// equilibration exists to honor. The scalar instantiation is the one the
// analytical case in this file already carries, so this adds none.
struct scalar_pose
{
    Eigen::Matrix<double, 1, 1> A;
    Eigen::Matrix<double, 1, 1> B;
};

auto common_scale_pose() -> scalar_pose
{
    scalar_pose pose;
    pose.A(0, 0) = 0.5;
    pose.B(0, 0) = 1.0;
    return pose;
}

auto weight_at(double scale) -> Eigen::Matrix<double, 1, 1>
{
    Eigen::Matrix<double, 1, 1> weight;
    weight(0, 0) = scale;
    return weight;
}

/// The gain the scalar pose implies, formed from the ratio of the returned
/// matrix to the common scale rather than at the caller's own scale. The gain is
/// homogeneous of degree zero, and the sum R + B'PB it is ordinarily read from
/// leaves the top of the range while the answer is still perfectly
/// representable, so forming that ratio first is what makes the oracle usable
/// across the whole range instead of only the middle of it.
auto scale_free_gain(const scalar_pose &pose, double P, double common_scale) -> double
{
    const double p = P / common_scale;
    return (pose.B(0, 0) * p * pose.A(0, 0)) / (1.0 + pose.B(0, 0) * p * pose.B(0, 0));
}

auto is_enumerated_dare_error(ctrlpp::dare_error error) -> bool
{
    switch(error)
    {
        case ctrlpp::dare_error::non_stabilizable:
        case ctrlpp::dare_error::non_finite_input:
        case ctrlpp::dare_error::singular_a:
        case ctrlpp::dare_error::singular_r:
        case ctrlpp::dare_error::singular_u11:
        case ctrlpp::dare_error::non_psd_solution:
        case ctrlpp::dare_error::schur_failed:
        case ctrlpp::dare_error::arithmetic_limit:
            return true;
    }
    return false;
}

struct sweep_tally
{
    std::size_t points{};
    std::size_t accepted{};
    std::size_t declined{};
    std::size_t compared{};
    std::size_t twin_declined{};
    std::size_t accepted_beyond_prior_reach{};
    std::size_t declined_beyond_ceiling{};
};

}

TEST_CASE("DARE carries a representable common-scaled pose", "[dare][hardening][precision]")
{
    // The pose that contradicted equilibration's own claim. Its unit-scale answer
    // is an ordinary normal double, and so is the answer at the common scale
    // below; the solve refused it because a squared quantity left the top of the
    // range while the evidence was being formed, not because the answer was
    // unrepresentable.
    const auto pose        = common_scale_pose();
    const auto unit_weight = weight_at(1.0);

    const auto unit = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, unit_weight, unit_weight);
    REQUIRE(unit.has_value());
    const auto unit_gain = ctrlpp::test::riccati_gain<double, 1, 1>(pose.A, pose.B, unit_weight, unit->P);

    const double relative_margin = ctrlpp::test::riccati_gain_agreement_margin<double, 1, 1>();

    const double common_scale = 1e155;
    const auto weight         = weight_at(common_scale);
    const auto scaled         = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight, weight);
    REQUIRE(scaled.has_value());

    // The answer is the unit-scale answer moved by the common factor.
    const double expected = unit->P(0, 0) * common_scale;
    CAPTURE(scaled->P(0, 0), expected);
    CHECK(std::abs(scaled->P(0, 0) - expected) <= relative_margin * expected);

    // The scale-invariant oracle. A residual is a necessary companion to this and
    // never a substitute: the gain is what two posings of the same problem have
    // in common as a mathematical identity.
    const auto scaled_gain = ctrlpp::test::riccati_gain<double, 1, 1>(pose.A, pose.B, weight, scaled->P);
    CAPTURE(unit_gain, scaled_gain);
    CHECK((unit_gain - scaled_gain).norm() <= relative_margin * unit_gain.norm());

    // Positive semi-definiteness holds at the scale actually returned, against
    // the shared pivot floor rather than a copy of it.
    Eigen::LDLT<Eigen::Matrix<double, 1, 1>> returned_psd(scaled->P);
    REQUIRE(returned_psd.info() == Eigen::Success);
    CHECK(returned_psd.vectorD().minCoeff() >= ctrlpp::detail::psd_pivot_floor<double, 1>(scaled->P));
}

TEST_CASE("DARE accepts inside the representable ceiling and declines outside it", "[dare][hardening][precision]")
{
    // Where the accepted range stops is a derived property, not a constant. The
    // returned matrix is the unit-scale answer times the common scale, so the
    // largest common scale whose answer is representable is the largest finite
    // value divided by the unit-scale answer's largest entry. Nothing here pins
    // the transition as an equality: acceptance is asserted inside the derived
    // ceiling and an enumerated decline outside it.
    const auto pose        = common_scale_pose();
    const auto unit_weight = weight_at(1.0);

    const auto unit = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, unit_weight, unit_weight);
    REQUIRE(unit.has_value());

    const double ceiling = std::numeric_limits<double>::max() / unit->P.cwiseAbs().maxCoeff();
    CAPTURE(ceiling);

    for(const double inside : {1e155, 1e308})
    {
        CAPTURE(inside);
        REQUIRE(inside < ceiling);
        const auto weight = weight_at(inside);
        const auto result = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight, weight);
        REQUIRE(result.has_value());
    }

    for(const double outside : {std::nextafter(ceiling, std::numeric_limits<double>::infinity()), std::numeric_limits<double>::max()})
    {
        CAPTURE(outside);
        REQUIRE(outside > ceiling);
        const auto weight = weight_at(outside);
        const auto result = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight, weight);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::dare_error::arithmetic_limit);
    }
}

TEST_CASE("DARE holds the common-scale identity across the representable range", "[dare][hardening][precision]")
{
    // The existing common-scale case stops an order of magnitude short of where
    // the defect lived, which is exactly why a green suite did not see it. This
    // one steps the common scale over every decade the type supports, in both
    // directions, at three fixed ratios, and holds every accepted point to the
    // gain the same pose implies at unit scale.
    //
    // The oracle is the gain of the pose the solver ACTUALLY RECEIVED, not of the
    // pose that was intended. Near the bottom of the range a weight formed as
    // `ratio * scale` can underflow, so the stored pair no longer realizes the
    // intended ratio; comparing it against the intended ratio's twin would report
    // a defect on a correct answer to the question actually asked. Dividing the
    // stored pair by its own input weighting is the same common rescale the
    // identity is about, so the twin is the right reference at every point,
    // faithful pose or not.
    const auto pose              = common_scale_pose();
    const double relative_margin = ctrlpp::test::riccati_gain_agreement_margin<double, 1, 1>();

    // The magnitude the previous verification stopped at on the equal-weight
    // pose, bisected to a relative width below 1e-13. Every accepted point above
    // it at that ratio is inside the region that implementation refused, which is
    // what makes the census non-vacuous. The census counts only that ratio,
    // because that is the ratio the number was measured at.
    const double prior_reach = 1.183617e154;

    sweep_tally tally;

    for(const double ratio : {1.0, 100.0, 0.01})
    {
        for(const bool ascending : {true, false})
        {
            for(int step = -323; step <= 308; ++step)
            {
                const int exponent        = ascending ? step : -step + (-323 + 308);
                const double common_scale = std::pow(10.0, static_cast<double>(exponent));
                CAPTURE(ratio, ascending, exponent, common_scale);

                const auto Q = weight_at(ratio * common_scale);
                const auto R = weight_at(common_scale);
                REQUIRE(R(0, 0) > 0.0);

                const auto result = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, Q, R);
                ++tally.points;
                if(!result)
                {
                    ++tally.declined;
                    CHECK(is_enumerated_dare_error(result.error()));
                    continue;
                }

                ++tally.accepted;
                if(ratio == 1.0 && common_scale > prior_reach)
                    ++tally.accepted_beyond_prior_reach;

                // The same pose divided by its own input weighting: one division
                // per entry, so the twin differs from an exact common rescale by
                // at most a single rounding, which the margin below absorbs.
                const auto twin = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight_at(Q(0, 0) / R(0, 0)), weight_at(1.0));
                if(!twin)
                {
                    ++tally.twin_declined;
                    continue;
                }

                ++tally.compared;
                const double reference_gain = scale_free_gain(pose, twin->P(0, 0), 1.0);
                const double implied_gain   = scale_free_gain(pose, result->P(0, 0), common_scale);
                CAPTURE(reference_gain, implied_gain, result->P(0, 0));
                CHECK(std::abs(implied_gain - reference_gain) <= relative_margin * std::abs(reference_gain));
            }
        }
    }

    // Every decade the type supports, at three ratios, in both directions.
    CHECK(tally.points == 3 * 632 * 2);
    CHECK(tally.accepted + tally.declined == tally.points);
    // Neither half of the census may pass vacuously: the sweep must have entered
    // the region the previous implementation refused, and it must have compared
    // real points rather than skipping them all.
    CHECK(tally.accepted_beyond_prior_reach > 0);
    CHECK(tally.compared > 0);
    CHECK(tally.twin_declined == 0);

    // The other half of the census: above the derived ceiling the answer is not
    // representable and the decline is enumerated. Nothing here pins where the
    // transition sits; it is a property of the type's range.
    const auto unit = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight_at(1.0), weight_at(1.0));
    REQUIRE(unit.has_value());
    const double ceiling = std::numeric_limits<double>::max() / unit->P.cwiseAbs().maxCoeff();
    for(const double outside : {std::nextafter(ceiling, std::numeric_limits<double>::infinity()), std::numeric_limits<double>::max()})
    {
        CAPTURE(outside);
        const auto weight = weight_at(outside);
        const auto result = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight, weight);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::dare_error::arithmetic_limit);
        ++tally.declined_beyond_ceiling;
    }
    CHECK(tally.declined_beyond_ceiling > 0);
}

TEST_CASE("DARE keeps the band whose only failing check is the gain's own scale", "[dare][hardening][precision]")
{
    // Near the top of the range the check at the CALLER's scale dissolves before
    // the answer does: the gain it recomputes needs `R + B'PB` formed there, and
    // that sum leaves the range while the solution is still an ordinary normal
    // number. That is `gain_unavailable` -- an absence of evidence with a known
    // cause -- and it is the one absence the solver accepts, because the
    // equilibrated verdict already covers the band and the answers in it are not
    // approximate. A rule that declined every unresolved verdict alike would
    // discard this band, and the census above would not notice: it only pins
    // that points ABOVE the derived ceiling decline.
    //
    // Nothing here is hardcoded. The ceiling is derived from the type's range,
    // the band's lower edge is found by walking down until the direct check
    // starts resolving again, and the correctness bar is exact equality against
    // the homogeneous truth rather than a tolerance.
    const auto pose = common_scale_pose();
    const auto unit = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight_at(1.0), weight_at(1.0));
    REQUIRE(unit.has_value());

    const double p_unit  = unit->P(0, 0);
    const double ceiling = std::numeric_limits<double>::max() / p_unit;

    int in_band = 0;
    for(double scale = ceiling; scale > 1.0; scale /= 2.0)
    {
        const auto weight = weight_at(scale);
        const auto result = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight, weight);
        REQUIRE(result.has_value());

        // Stop at the lower edge: below it the direct check resolves again and
        // the pose is no longer evidence about this rule.
        if(ctrlpp::detail::verify_dare_solution<double, 1, 1>(pose.A, pose.B, weight, weight, result->P) != ctrlpp::detail::dare_verification::gain_unavailable)
            break;

        ++in_band;
        // Homogeneity of degree one: the answer is exactly `scale * p_unit`, and
        // the rescale is a single multiplication, so the returned matrix carries
        // the unit answer's every bit. Checked by division to stay in range.
        CAPTURE(scale, result->P(0, 0), p_unit);
        CHECK(result->P(0, 0) / scale == p_unit);
    }

    // The band exists and the walk entered it, so neither half passes vacuously.
    CHECK(in_band > 0);
}

TEST_CASE("DARE original-scale residual guard still asserts where a squared magnitude would not", "[dare][hardening][precision]")
{
    // The guard integrity half of the change, and it is independent of reach. At
    // a common scale of 1e-200 every term of the residual expression has a sum of
    // squares below the smallest normal value, so a plain Frobenius magnitude
    // gives a residual of zero against a scale of zero and the acceptance test
    // becomes `0 <= 0` -- an unconditionally passing guard over roughly a hundred
    // and twenty decades. Resolving both magnitudes through the largest-entry
    // rescale keeps them meaningful, which is checked here by handing the
    // verification a matrix that does NOT solve the equation and requiring it to
    // say so at exactly that scale.
    const auto pose           = common_scale_pose();
    const double common_scale = 1e-200;
    const auto weight         = weight_at(common_scale);

    const auto solved = ctrlpp::dare<double, 1, 1>(pose.A, pose.B, weight, weight);
    REQUIRE(solved.has_value());

    // The correct answer verifies at its own scale.
    CHECK(ctrlpp::detail::verify_dare_solution<double, 1, 1>(pose.A, pose.B, weight, weight, solved->P) == ctrlpp::detail::dare_verification::verified);

    // A positive semi-definite matrix of the same magnitude that solves nothing
    // is refuted. Under the plain magnitude both sides of the comparison were
    // exactly zero and this matrix passed.
    auto wrong = solved->P;
    wrong(0, 0) *= 2.0;
    CHECK(ctrlpp::detail::verify_dare_solution<double, 1, 1>(pose.A, pose.B, weight, weight, wrong) == ctrlpp::detail::dare_verification::refuted);

    // And the degenerate reading is gone at the source: the residual scale this
    // pose produces is nonzero, where a plain sum of squares gives exactly zero.
    const auto plain_scale    = solved->P.norm();
    const auto resolved_scale = ctrlpp::detail::resolve_magnitude(solved->P);
    CAPTURE(plain_scale, resolved_scale.value);
    CHECK(plain_scale == 0.0);
    REQUIRE(resolved_scale.resolved);
    CHECK(resolved_scale.value > 0.0);
}

TEST_CASE("DARE keeps the two poses its randomized oracle used to abort on", "[dare][hardening][precision]")
{
    // Both of these were reached by the discrete randomized target and both made
    // it abort, at 1.078 and 1.124 times the tolerance that oracle then carried.
    // Both answers are RIGHT: against an independent extended-precision solution
    // of the identical pose, computed by Newton-Kleinman rather than by a Schur
    // decomposition, they agree to more than nine decimal digits -- comfortably
    // more than half of binary64's significand. The oracle was aborting on
    // correct answers, and its verdict at 1.078x was inside the reproducibility
    // noise of the quantity it compared: the identical source built with fused
    // multiply-add contraction returns a MORE accurate answer on both poses and
    // does not abort at all.
    //
    // They are kept here as correct-answer regressions, with their exact entries,
    // so the band between two disagreeing bounds cannot silently reopen. Note
    // what is asserted: not a residual, which cannot tell a right answer from a
    // wrong one on this population, but the retained accuracy of the answer.
    auto keeps_more_than_half_the_significand =
        [](const Eigen::Matrix<double, 2, 2>& A, const Eigen::Matrix<double, 2, 1>& B,
           const Eigen::Matrix<double, 2, 2>& Q, const Eigen::Matrix<double, 1, 1>& R)
    {
        const auto solved = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
        REQUIRE(solved.has_value());

        // The one rule, quoted. Both poses keep more than half the significand,
        // so both are accepted, and the shared verdict says so directly rather
        // than by inference from the solver having returned a value.
        CHECK(ctrlpp::test::riccati_accuracy_of<double, 2, 1>(A, B, Q, R, solved->P)
              == ctrlpp::detail::riccati_accuracy::within);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(solved->P);
        for(int i = 0; i < 2; ++i)
            CHECK(pes.eigenvalues()(i) > 0.0);
    };

    Eigen::Matrix<double, 2, 2> first_A;
    first_A << -0.4509803921568627, -0.12401960784313724, 2.0, -2.0;
    Eigen::Matrix<double, 2, 1> first_B;
    first_B << -2.0, -2.0;
    Eigen::Matrix<double, 2, 2> first_Q;
    first_Q << 8.0, 4.0, 4.0, 4.0;
    Eigen::Matrix<double, 1, 1> first_R;
    first_R << 0.004;
    keeps_more_than_half_the_significand(first_A, first_B, first_Q, first_R);

    Eigen::Matrix<double, 2, 2> second_A;
    second_A << 0.0, 0.42742921411995211, 2.0, -2.0;
    Eigen::Matrix<double, 2, 1> second_B;
    second_B << -2.0, -1.999999231413548;
    Eigen::Matrix<double, 2, 2> second_Q;
    second_Q << 0.030762027265631202, -2.3674744784873969e-06, -2.3674744784873969e-06, 8.0;
    Eigen::Matrix<double, 1, 1> second_R;
    second_R << 0.004;
    keeps_more_than_half_the_significand(second_A, second_B, second_Q, second_R);
}

TEST_CASE("DARE refuses an answer that has lost more than half the significand", "[dare][hardening][precision]")
{
    // This pose is outside every filter the randomized target applies, and it is
    // the failure the residual margin could not see. The solver's own answer here
    // is wrong in the twentieth bit: against an independent extended-precision
    // solution by a different algorithm its relative forward error is 5.7e-07,
    // which is thirty-eight times past the half-significand line, while the
    // residual that answer produces sits comfortably INSIDE the counted-operation
    // envelope the postcondition used to compare against. A residual bound
    // therefore certified it, and no threshold on the residual could have done
    // otherwise: on this population the answers that keep half the significand
    // and the answers that lose it are contiguous on that quantity.
    //
    // The verdict is stable under instruction selection, which is why this pose
    // and not a marginal one is the regression: it stays past the line under six
    // optimization settings across two compilers, including the fused-multiply-add
    // build that flips other candidates.
    Eigen::Matrix<double, 2, 2> A;
    A << -1.3314095436239524, -1.3025433868459988, -2.0, -2.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 2.0, -2.0;
    Eigen::Matrix<double, 2, 2> Q;
    Q << 0.12179221216883218, -0.48675460347327504, -0.48675460347327504, 4.2098620441168375;
    Eigen::Matrix<double, 1, 1> R;
    R << 0.05248182488128475;

    const auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::arithmetic_limit);
}
