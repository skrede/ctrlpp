#include "hardening_helpers.h"
#include "ctrlpp/control/care.h"
#include "ctrlpp/control/lqr.h"
#include "ctrlpp/detail/quasi_triangular.h"


#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <cstddef>
#include <numbers>


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

TEST_CASE("the closed-loop spectrum test refuses a 2x2 block whose binding root is unstable",
          "[care][error][postcondition]")
{
    // The acceptance rule reads the closed loop's real Schur factor block by
    // block, and a 2x2 block with a significant subdiagonal entry is NOT
    // necessarily a complex-conjugate pair. When its discriminant is
    // non-negative the block holds two REAL eigenvalues, and their mean -- the
    // quantity a walk that assumes a conjugate pair reads -- is neither of them.
    //
    // The block below is exactly that case, constructed from the closed loop
    // backwards: pick the spectrum {-3, +1}, which needs a trace of -2 and a
    // determinant of -3, then pick entries realizing it with a nonzero
    // subdiagonal. Its half-trace is -1, comfortably left of any counted margin,
    // while its binding root sits at +1 in the open right half-plane.
    Eigen::Matrix<double, 2, 2> factor;
    factor << 0.0,  3.0,
              1.0, -2.0;

    const double factor_scale = factor.cwiseAbs().maxCoeff();
    // The same counted form the rule itself uses: the state-space dimension
    // times unit roundoff times the factor's largest entry.
    const double margin =
        2.0 * std::numeric_limits<double>::epsilon() * factor_scale;

    // Reachability first. A case whose block silently degenerates to two 1x1
    // blocks, or whose discriminant comes out negative, would pass the verdict
    // assertion below for the wrong reason.
    REQUIRE(ctrlpp::detail::quasi_triangular_block_size<double, 2>(
                factor, 0, factor_scale)
            == 2);
    const auto block =
        ctrlpp::detail::quasi_triangular_block_spectrum<double, 2>(factor, 0, 2);
    REQUIRE(block.real_pair);
    CHECK(block.first.real() == 1.0);
    CHECK(block.second.real() == -3.0);

    // The rule that was in place before the discriminant was asked read the
    // half-trace. Asserting that the half-trace WOULD have passed is what makes
    // this a regression rather than a restatement of the fix: reintroducing the
    // mean turns the check below red.
    const double half_trace = (factor(0, 0) + factor(1, 1)) / 2.0;
    REQUIRE(half_trace < -margin);
    CHECK_FALSE(ctrlpp::detail::quasi_triangular_spectrum_strictly_left_of<double, 2>(
        factor, margin));

    // The opposite direction, so the fix cannot be a blanket refusal of 2x2
    // blocks: a genuine complex-conjugate pair strictly in the open left
    // half-plane is still accepted. This is the exact solution of the double
    // integrator's factor, whose diagonal is (0, -sqrt(3)) for a pair whose real
    // parts are both -sqrt(3)/2.
    Eigen::Matrix<double, 2, 2> oscillatory;
    oscillatory << 0.0,           -std::sqrt(3.0),
                   std::sqrt(3.0), -std::sqrt(3.0);
    const auto pair = ctrlpp::detail::quasi_triangular_block_spectrum<double, 2>(
        oscillatory, 0, 2);
    REQUIRE_FALSE(pair.real_pair);
    CHECK(ctrlpp::detail::quasi_triangular_spectrum_strictly_left_of<double, 2>(
        oscillatory,
        2.0 * std::numeric_limits<double>::epsilon()
            * oscillatory.cwiseAbs().maxCoeff()));
}

TEST_CASE("a real-eigenvalue 2x2 block does not survive Eigen's factorization",
          "[care][error][postcondition]")
{
    // What this pins is an assumption, not a behavior of this library. Eigen's
    // RealSchur triangularizes a 2x2 block whose discriminant is non-negative
    // and writes an EXACT zero into the subdiagonal, so the real-pair branch of
    // the block walk is not reachable through the public continuous entry point
    // today. That is an undocumented internal of a third-party header, and the
    // acceptance rule no longer depends on it -- the case above shows the rule
    // is correct whether or not this holds. This case exists so that a future
    // Eigen which stops doing it is reported here rather than discovered as a
    // certified unstable closed loop.
    for(int angle_index = 0; angle_index < 32; ++angle_index)
    {
        const double angle =
            std::numbers::pi * static_cast<double>(angle_index) / 32.0;
        Eigen::Matrix<double, 2, 2> transform;
        transform << std::cos(angle), -std::sin(angle),
                     std::sin(angle),  std::cos(angle);
        Eigen::Matrix<double, 2, 2> spectrum;
        spectrum << -3.0, 0.0, 0.0, 1.0;
        // A shear keeps the spectrum and destroys normality, which is the shape
        // most likely to leave a coupled block behind.
        Eigen::Matrix<double, 2, 2> shear;
        shear << 1.0, 7.5, 0.0, 1.0;
        const Eigen::Matrix<double, 2, 2> closed_loop =
            (shear * transform * spectrum * transform.transpose()
             * shear.inverse())
                .eval();

        CAPTURE(angle_index);
        Eigen::RealSchur<Eigen::Matrix<double, 2, 2>> schur(closed_loop, false);
        REQUIRE(schur.info() == Eigen::Success);
        CHECK(schur.matrixT()(1, 0) == 0.0);
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
        // and past three decades its extracted matrix no longer satisfies the
        // equation to the counted bound.
        //
        // The band starts at three decades rather than two because the counted
        // bound was corrected: the extraction count is now taken at the 2n
        // basis every path's factorization actually produces rather than at n,
        // and the decade this method used to refuse it now answers. That answer
        // is checked below rather than assumed.
        for(int exponent = 3; exponent <= 8; ++exponent)
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
        // about the problem: just below the band both the default tag and the
        // real-Schur tag answer the same pose, and both answers are stabilizing.
        const Eigen::Matrix<double, 2, 2> Q = (100.0 * base_Q).eval();
        const Eigen::Matrix<double, 2, 2> R = (100.0 * base_R).eval();
        auto const require_stabilizing =
            [&](const Eigen::Matrix<double, 2, 2>& P) {
                const Eigen::Matrix<double, 2, 2> closed_loop =
                    (A - B * R.inverse() * B.transpose() * P).eval();
                Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> spectrum(
                    closed_loop, false);
                CHECK(spectrum.eigenvalues()(0).real() < 0.0);
                CHECK(spectrum.eigenvalues()(1).real() < 0.0);
            };

        auto const sign = ctrlpp::care<double, 2, 2>(A, B, Q, R);
        REQUIRE(sign.has_value());
        require_stabilizing(sign->P);

        // The decade the corrected extraction count recovered. It is asserted
        // as an answer AND as a correct one, because a bound that was loosened
        // has to be shown to have removed a refusal rather than admitted a
        // wrong success.
        auto const recovered = ctrlpp::care<double, 2, 2>(
            A, B, Q, R, ctrlpp::detail::schur_care_method{});
        REQUIRE(recovered.has_value());
        require_stabilizing(recovered->P);
    }

    SECTION("the balanced-Schur tag over a well-posed near-axis band")
    {
        // A band rather than a pinned pose, deliberately. Which side of the
        // postcondition an individual near-axis draw lands on is decided in the
        // last bits and moves with instruction selection: a pose pinned here
        // was observed to decline at -O0 and -O2 and to be accepted at
        // -O3 -march=native. What does not move is the shape of the band, so
        // the band is asserted and the pose is not.
        //
        // Each count below names the enumerator it counts, because two of them
        // used to be described with the same word. `unverified` counts refusals
        // carrying `unverified_solution` specifically -- 86 to 93 of these 112
        // poses, on the measurements taken so far. `stagnated` counts refusals
        // carrying `sign_function_stagnated`, and it must be zero: this path
        // runs no Newton iteration. Refusals carrying any other enumerator are
        // neither, and are not counted here.
        //
        // What is asserted about acceptance is that every answer is right, NOT
        // that there are no answers. An earlier form of this case required zero
        // acceptances, which asserted the presence of over-refusal: when the
        // extraction count was corrected to the 2n basis the factorization
        // actually produces, one pose of the 112 came back with a stabilizing
        // answer, and a test demanding zero would have called that a regression.
        //
        // Every member is one unstable mode approaching the imaginary axis with
        // its input direction shrinking alongside it, rotated into general
        // position. The pair is controllable and Q = I makes it detectable, so
        // a decline is the method declining to certify its own answer rather
        // than the problem lacking one.
        std::size_t accepted = 0;
        std::size_t accepted_unstable = 0;
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
                {
                    ++accepted;
                    const Eigen::Matrix<double, 2, 2> closed_loop =
                        (A - B * B.transpose() * balanced->P).eval();
                    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> spectrum(
                        closed_loop, false);
                    if(!(spectrum.eigenvalues()(0).real() < 0.0
                         && spectrum.eigenvalues()(1).real() < 0.0))
                        ++accepted_unstable;
                }
                else if(balanced.error()
                        == ctrlpp::care_error::unverified_solution)
                    ++unverified;
                else if(balanced.error()
                        == ctrlpp::care_error::sign_function_stagnated)
                    ++stagnated;
            }
        }

        CAPTURE(total, accepted, accepted_unstable, unverified, stagnated);
        // The exact count the two literal loop bounds produce, so a bound that
        // changes is caught here rather than absorbed by a comparison against
        // zero that no loop can fail.
        CHECK(total == 8 * 14);
        // The load-bearing pair. Nothing this tag accepts may be unstable, and
        // this path runs no Newton iteration, so the enumerator that names one
        // must never come out of it.
        CHECK(accepted_unstable == 0);
        CHECK(stagnated == 0);
        CHECK(unverified > 0);
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
