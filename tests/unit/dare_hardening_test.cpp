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

TEST_CASE("DARE accepts an R that is ill-conditioned but not singular", "[dare][hardening][robustness]")
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

    constexpr double eps              = std::numeric_limits<double>::epsilon();
    const double gain_relative_margin = std::sqrt(ctrlpp::test::riccati_residual_ops<2, 1> * eps);
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
    constexpr double eps              = std::numeric_limits<double>::epsilon();
    const double gain_relative_margin = std::sqrt(ctrlpp::test::riccati_residual_ops<2, 1> * eps);
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

    constexpr double eps         = std::numeric_limits<double>::epsilon();
    const double relative_margin = std::sqrt(ctrlpp::test::riccati_residual_ops<1, 1> * eps);

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
    constexpr double eps         = std::numeric_limits<double>::epsilon();
    const double relative_margin = std::sqrt(ctrlpp::test::riccati_residual_ops<1, 1> * eps);

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
