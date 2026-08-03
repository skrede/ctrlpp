#include "hardening_helpers.h"
#include "ctrlpp/control/care.h"
#include "ctrlpp/control/lqr.h"


#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>


TEST_CASE("CARE non-LHP-stabilizable system fails with non_lhp_stabilizable or singular_u11",
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
    CHECK((result.error() == ctrlpp::care_error::non_lhp_stabilizable
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

TEST_CASE("CARE cross-weight overload classifies every non-finite operand",
          "[care][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0, -0.5, -0.3;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Matrix<double, 2, 1> N;
    N << 0.1, 0.2;

    auto require_non_finite = [](auto const& a, auto const& b, auto const& q,
                                 auto const& r, auto const& n) {
        auto const result = ctrlpp::care<double, 2, 1>(a, b, q, r, n);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::care_error::non_finite_input);
    };

    auto bad_A = A;
    bad_A(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require_non_finite(bad_A, B, Q, R, N);

    auto bad_B = B;
    bad_B(0, 0) = std::numeric_limits<double>::infinity();
    require_non_finite(A, bad_B, Q, R, N);

    auto bad_Q = Q;
    bad_Q(0, 0) = -std::numeric_limits<double>::infinity();
    require_non_finite(A, B, bad_Q, R, N);

    auto bad_R = R;
    bad_R(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require_non_finite(A, B, Q, bad_R, N);

    auto bad_N = N;
    bad_N(0, 0) = std::numeric_limits<double>::infinity();
    require_non_finite(A, B, Q, R, bad_N);
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
        || result.error() == ctrlpp::care_error::non_lhp_stabilizable
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
            || result.error() == ctrlpp::care_error::non_lhp_stabilizable
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

TEST_CASE("CARE refuses a singular R rather than naming a symptom of it",
          "[care][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0, 0.0, 0.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();

    SECTION("a zero weighting")
    {
        Eigen::Matrix<double, 1, 1> R;
        R << 0.0;
        REQUIRE(R.allFinite());

        // The Hamiltonian build needs R^{-1} for B R^{-1} B^T. Before the
        // enumerator existed this reported non_finite_input, a false statement
        // about data the caller chose.
        auto const result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::care_error::singular_r);

        // The Schur path assembles the same Hamiltonian, so it refuses identically.
        auto const schur = ctrlpp::care<double, 2, 1>(A, B, Q, R, ctrlpp::detail::schur_care_method{});
        REQUIRE_FALSE(schur.has_value());
        CHECK(schur.error() == ctrlpp::care_error::singular_r);

        // The cross-weight overload inverts R before the build sees it and
        // carries its own copy of the test.
        Eigen::Matrix<double, 2, 1> N;
        N << 0.1, 0.2;
        auto const crossed = ctrlpp::care<double, 2, 1>(A, B, Q, R, N);
        REQUIRE_FALSE(crossed.has_value());
        CHECK(crossed.error() == ctrlpp::care_error::singular_r);
    }

    SECTION("a rank-deficient but nonzero weighting")
    {
        // The dangerous half: this one never goes non-finite. The rank-revealing
        // QR solve returns a least-squares answer over the leading rank columns,
        // so without an explicit rank test the solver would report SUCCESS on a
        // Hamiltonian that is not the one the problem defines.
        Eigen::Matrix<double, 2, 2> B2;
        B2 << 0.0, 0.0, 1.0, 1.0;
        auto R2 = Eigen::Matrix<double, 2, 2>::Zero().eval();
        R2(0, 0) = 1.0;
        REQUIRE(R2.allFinite());

        auto const result = ctrlpp::care<double, 2, 2>(A, B2, Q, R2);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::care_error::singular_r);
    }
}

TEST_CASE("CARE reports an unverifiable extracted result with a method-neutral cause",
          "[care][error]")
{
    // Both Schur variants used to return whatever their extraction produced.
    // They now hold it to the same postconditions the default path does, and a
    // failure of those postconditions is reported as `unverified_solution` --
    // never as `sign_function_stagnated`, whose four documented cases are all
    // statements about a Newton iteration that these paths do not run.
    SECTION("the real-Schur tag on a common weight rescale it cannot certify")
    {
        Eigen::Matrix<double, 2, 2> A;
        A << -1.1, 0.3, -0.2, -1.4;
        Eigen::Matrix<double, 2, 2> B;
        B << 0.8, 0.1, -0.15, 0.65;
        Eigen::Matrix<double, 2, 2> base_Q;
        base_Q << 1.0, 0.2, 0.2, 1.7;
        Eigen::Matrix<double, 2, 2> base_R;
        base_R << 1.3, 0.1, 0.1, 0.9;

        // A common rescale of Q and R leaves the gain and the closed-loop
        // spectrum where they were, so every pose in this band has the same
        // answer up to the scale. What the rescale does move is the
        // conditioning of the invariant subspace this method extracts from,
        // and past two decades its extracted matrix no longer satisfies the
        // equation to the counted bound.
        for(int exponent = 2; exponent <= 8; ++exponent)
        {
            const double scale = std::pow(10.0, static_cast<double>(exponent));
            const Eigen::Matrix<double, 2, 2> Q = (scale * base_Q).eval();
            const Eigen::Matrix<double, 2, 2> R = (scale * base_R).eval();

            CAPTURE(exponent);
            auto const schur = ctrlpp::care<double, 2, 2>(
                A, B, Q, R, ctrlpp::detail::schur_care_method{});
            REQUIRE_FALSE(schur.has_value());
            CHECK(schur.error() == ctrlpp::care_error::unverified_solution);
            CHECK(schur.error() != ctrlpp::care_error::sign_function_stagnated);
        }

        // The refusal is a statement about that method's extraction and not
        // about the problem: at the bottom of the band the default tag answers
        // the same pose, and its answer is stabilizing.
        const Eigen::Matrix<double, 2, 2> Q = (100.0 * base_Q).eval();
        const Eigen::Matrix<double, 2, 2> R = (100.0 * base_R).eval();
        auto const sign = ctrlpp::care<double, 2, 2>(A, B, Q, R);
        REQUIRE(sign.has_value());
        const Eigen::Matrix<double, 2, 2> closed_loop =
            (A - B * R.inverse() * B.transpose() * sign->P).eval();
        Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> spectrum(closed_loop,
                                                                 false);
        CHECK(spectrum.eigenvalues()(0).real() < 0.0);
        CHECK(spectrum.eigenvalues()(1).real() < 0.0);
    }

    SECTION("the balanced-Schur tag over a well-posed near-axis band")
    {
        // A band rather than a pinned pose, deliberately. Which side of the
        // postcondition an individual near-axis draw lands on is decided in the
        // last bits and moves with instruction selection: a pose pinned here
        // was observed to decline at -O0 and -O2 and to be accepted at
        // -O3 -march=native. What does not move is the density -- 86 to 93 of
        // these 112 poses decline, none is accepted, and none reports the
        // sign-function cause -- so the band is asserted and the pose is not.
        //
        // Every member is one unstable mode approaching the imaginary axis with
        // its input direction shrinking alongside it, rotated into general
        // position. The pair is controllable and Q = I makes it detectable, so
        // a decline is the method declining to certify its own answer rather
        // than the problem lacking one.
        std::size_t accepted = 0;
        std::size_t unverified = 0;
        std::size_t stagnated = 0;
        std::size_t total = 0;

        for(int angle_index = 1; angle_index <= 8; ++angle_index)
        {
            const double angle = static_cast<double>(angle_index) / 8.0;
            const double cosine = std::cos(angle);
            const double sine = std::sin(angle);
            Eigen::Matrix<double, 2, 2> transform;
            transform << cosine, -sine, sine, cosine;

            for(int exponent = 2; exponent <= 15; ++exponent)
            {
                const double delta =
                    std::pow(10.0, -static_cast<double>(exponent));
                Eigen::Matrix<double, 2, 2> diagonal_a;
                diagonal_a << 0.5 * delta, 0.0, 0.0, -0.8;
                Eigen::Matrix<double, 2, 2> diagonal_b;
                diagonal_b << std::sqrt(0.75) * delta, 0.0, 0.0, 0.6;
                const Eigen::Matrix<double, 2, 2> A =
                    (transform * diagonal_a * transform.transpose()).eval();
                const Eigen::Matrix<double, 2, 2> B =
                    (transform * diagonal_b).eval();
                auto const I = Eigen::Matrix<double, 2, 2>::Identity();

                auto const balanced = ctrlpp::care<double, 2, 2>(
                    A, B, I, I,
                    ctrlpp::detail::balanced_schur_care_method{});
                ++total;
                if(balanced.has_value())
                    ++accepted;
                else if(balanced.error()
                        == ctrlpp::care_error::unverified_solution)
                    ++unverified;
                else if(balanced.error()
                        == ctrlpp::care_error::sign_function_stagnated)
                    ++stagnated;
            }
        }

        CAPTURE(total, accepted, unverified, stagnated);
        CHECK(total > 0);
        CHECK(accepted == 0);
        CHECK(unverified > 0);
        // The load-bearing one: this path runs no Newton iteration, so the
        // enumerator that names one must never come out of it.
        CHECK(stagnated == 0);
    }
}

TEST_CASE("lqr_gain_continuous refuses a singular R its own factorization hid",
          "[lqr][continuous][error]")
{
    // This surface forms R^{-1} through LDLT rather than through the Hamiltonian
    // build, and that factorization fails QUIETLY: its solve zeroes the
    // rank-deficient directions instead of producing infinities, so a singular R
    // yielded a finite R^{-1} of zeros, a finite Hamiltonian describing a plant
    // with no control authority, and a sign-function iteration that stagnated on
    // it. Nothing on that path was ever non-finite, so no downstream check could
    // catch it and the observed enumerator was sign_function_stagnated -- which
    // sends the caller to look at convergence.
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0, 0.0, 0.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.0;

    auto const result = ctrlpp::lqr_gain_continuous<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::care_error::singular_r);

    // The boundary the test must not overshoot: a small but nonsingular
    // weighting is a conditioning question, not a domain violation.
    Eigen::Matrix<double, 1, 1> R_small;
    R_small << 1e-9;
    auto const accepted = ctrlpp::lqr_gain_continuous<double, 2, 1>(A, B, Q, R_small);
    REQUIRE(accepted.has_value());
}
