// What the oracles in this file decide.
//
// These cases pin the rotation primitives to the answers their own algebra
// fixes, and nothing weaker. Three shapes appear, and which one a case uses is
// a statement about the mathematics, not a stylistic choice:
//
//  * Exact equality, where the answer is representable AND the computation
//    necessarily reaches it. The logarithm of a quaternion whose vector part is
//    exactly zero is the exact zero vector, because the Taylor branch returns a
//    scalar multiple of that vector; no tolerance can make that statement
//    stronger and any tolerance makes it weaker.
//  * A counted-operation budget, where a rounding genuinely enters. Every
//    budget below is a named integer constant whose count is enumerated in the
//    comment above it, multiplied by the scalar type's epsilon and by the
//    largest operand that entered. No bare numeric tolerance appears anywhere
//    in this file.
//  * A typed rejection, for the two inputs that have no unit representative.
//
// What they deliberately do not decide. Nothing here fixes the accuracy of the
// underlying transcendental functions beyond the counted budgets, and nothing
// here asserts a hemisphere convention for the logarithm beyond the
// canonicalization the function documents. The exp/log round trip is exercised
// on both sides of the small-angle branch threshold so the branch choice is
// covered, but the two branches agree to the last bit at an angle of 1e-12, so
// a case at that scale pins the value and NOT the branch -- which is why the
// threshold sweep exists as a separate case.

#include "hardening_helpers.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double so3_eps = std::numeric_limits<double>::epsilon();

// Rounded operations along the longest chain from a rotation vector to one
// entry of the rotation matrix it generates. Enumerated: the exponential forms
// the angle as a three-term norm (three multiplies, two additions, one square
// root = six), halves it (one), takes one sine and one cosine (two), divides
// the sine by the angle (one) and scales one vector component (one) -- ten to a
// quaternion coefficient. The quaternion-to-matrix conversion then forms one
// entry from two products of doubled coefficients and up to two sums (four).
constexpr int so3_exp_to_matrix_ops = 14;

// The same chain, but ending at the quaternion's own norm rather than at a
// matrix entry: ten to each coefficient, then four squares, three sums and one
// square root to combine them.
constexpr int so3_exp_norm_ops = 18;

// Rounded operations along an exp followed by a log. The ten above, then the
// logarithm's vector norm (three multiplies, two additions, one square root =
// six), one arc tangent, one division by the vector norm, one doubling and one
// scaling of the component -- sixteen in total.
constexpr int so3_exp_log_roundtrip_ops = 16;

}

// ── SO(3) hardening ────────────────────────────────────────────────────────────

TEST_CASE("SO3 zero quaternion", "[so3][hardening][negative]")
{
    Eigen::Quaterniond q;
    q.w() = 0.0;
    q.x() = 0.0;
    q.y() = 0.0;
    q.z() = 0.0;

    // The logarithm of the zero quaternion is the EXACT zero vector, and that is
    // a fact about the branch structure rather than a coincidence: the vector
    // part has norm zero, which selects the Taylor branch, whose scale factor is
    // the constant two, and two times the exact zero vector is the exact zero
    // vector. An implementation that took the arc-tangent branch here would
    // divide by zero and return NaN, so this is falsifiable and exact.
    auto v = ctrlpp::so3::log(q);
    REQUIRE(v == ctrlpp::Vector<double, 3>::Zero());

    // The zero quaternion carries no direction, so no scaling of it lands on
    // the unit sphere. normalize reports that rather than handing the input
    // back: the linear-algebra library's normalizing member returns a COPY of
    // its argument when squaredNorm() > 0 does not hold (Eigen 3.4.0,
    // Eigen/src/Core/Dot.h:122-134), which is what a unit-norm postcondition
    // built on it would silently fail to deliver here.
    auto const qn = ctrlpp::so3::normalize(q);
    REQUIRE_FALSE(qn.has_value());
    REQUIRE(qn.error() == ctrlpp::so3_error::zero_quaternion);
}

TEST_CASE("SO3 normalize delivers unit norm or names why it cannot", "[so3][hardening][negative]")
{
    using ctrlpp::so3_error;
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();

    // A counted-operation budget for the realized norm: the normalization
    // performs one division by the largest coefficient, one four-term sum of
    // squares, one square root, and one division by that root, and the norm
    // measured here adds a second sum of squares and square root -- six
    // rounded operations on the coefficient that carries the result.
    constexpr int normalize_norm_ops = 6;
    constexpr double norm_budget = normalize_norm_ops * std::numeric_limits<double>::epsilon();

    SECTION("a non-finite coefficient is refused, not propagated")
    {
        for(auto const bad : {nan, inf, -inf})
        {
            auto const r = ctrlpp::so3::normalize(Eigen::Quaterniond{bad, 0.0, 0.0, 0.0});
            REQUIRE_FALSE(r.has_value());
            REQUIRE(r.error() == so3_error::non_finite_input);

            auto const r2 = ctrlpp::so3::normalize(Eigen::Quaterniond{1.0, bad, 0.0, 0.0});
            REQUIRE_FALSE(r2.has_value());
            REQUIRE(r2.error() == so3_error::non_finite_input);
        }
    }

    SECTION("an ordinary non-unit quaternion is normalized")
    {
        auto const r = ctrlpp::so3::normalize(Eigen::Quaterniond{2.0, 0.0, 0.0, 0.0});
        REQUIRE(r.has_value());
        REQUIRE(r->w() == 1.0);
        REQUIRE_THAT(r->norm(), WithinAbs(1.0, norm_budget));
    }

    SECTION("a quaternion whose squared norm underflows is normalized, not returned unchanged")
    {
        // 1e-200 squares to zero in double, so squaredNorm() > 0 fails on a
        // finite quaternion with a perfectly well-defined direction.
        auto const q = Eigen::Quaterniond{1e-200, 2e-200, 0.0, 0.0};
        REQUIRE(q.coeffs().squaredNorm() == 0.0);

        auto const r = ctrlpp::so3::normalize(q);
        REQUIRE(r.has_value());
        REQUIRE_THAT(r->norm(), WithinAbs(1.0, norm_budget));
        REQUIRE_THAT(r->w() / r->x(), WithinAbs(0.5, norm_budget));
    }

    SECTION("a quaternion whose squared norm overflows is normalized, not collapsed to zero")
    {
        // 1e200 squares to infinity, so the divisor is infinite and a division
        // formed from it yields the zero quaternion.
        auto const q = Eigen::Quaterniond{1e200, 2e200, 0.0, 0.0};
        REQUIRE(std::isinf(q.coeffs().squaredNorm()));

        auto const r = ctrlpp::so3::normalize(q);
        REQUIRE(r.has_value());
        REQUIRE_THAT(r->norm(), WithinAbs(1.0, norm_budget));
        REQUIRE_THAT(r->w() / r->x(), WithinAbs(0.5, norm_budget));
    }
}

TEST_CASE("SO3 non-unit quaternion", "[so3][hardening][negative]")
{
    Eigen::Quaterniond q;
    q.w() = 2.0;
    q.x() = 0.0;
    q.y() = 0.0;
    q.z() = 0.0;

    // The scalar-only quaternion has an exactly zero vector part, so the
    // logarithm takes the Taylor branch and returns twice that vector: the exact
    // zero rotation. Finiteness was true of infinitely many wrong answers; this
    // is the only right one, and it says what "non-unit" costs the caller here,
    // namely nothing at all, because the logarithm reads only the direction.
    auto v = ctrlpp::so3::log(q);
    REQUIRE(v == ctrlpp::Vector<double, 3>::Zero());
}

TEST_CASE("SO3 90-degree rotation around z-axis", "[so3][hardening][precision]")
{
    ctrlpp::Vector<double, 3> phi;
    phi << 0.0, 0.0, std::numbers::pi / 2.0;

    auto q = ctrlpp::so3::exp(phi);
    auto R = q.toRotationMatrix();

    // Every entry is an order-one quantity, so the operand scale is one.
    constexpr double matrix_budget = so3_exp_to_matrix_ops * so3_eps;

    // Expected: [[0,-1,0],[1,0,0],[0,0,1]]
    REQUIRE_THAT(R(0, 0), WithinAbs(0.0, matrix_budget));
    REQUIRE_THAT(R(0, 1), WithinAbs(-1.0, matrix_budget));
    REQUIRE_THAT(R(1, 0), WithinAbs(1.0, matrix_budget));
    REQUIRE_THAT(R(1, 1), WithinAbs(0.0, matrix_budget));
    REQUIRE_THAT(R(2, 2), WithinAbs(1.0, matrix_budget));
}

TEST_CASE("SO3 180-degree rotation (near singularity)", "[so3][hardening][robustness]")
{
    ctrlpp::Vector<double, 3> phi;
    phi << std::numbers::pi, 0.0, 0.0;

    auto q = ctrlpp::so3::exp(phi);
    REQUIRE_THAT(q.norm(), WithinAbs(1.0, so3_exp_norm_ops * so3_eps));

    // Half a turn is the branch cut of the logarithm: the scalar part is within
    // an ulp of zero there, which is exactly where the arc tangent's own
    // conditioning is worst. The budget is therefore relative to pi, the operand
    // that entered, and counted along the round trip rather than guessed.
    auto recovered = ctrlpp::so3::log(q);
    REQUIRE_THAT(recovered.norm(),
                 WithinRel(std::numbers::pi, so3_exp_log_roundtrip_ops * so3_eps));

    // The off-axis components are exact: the exponential scaled two exactly zero
    // vector components, and the logarithm scaled them again, so no rounding can
    // move them off zero. A logarithm that leaked the rotation axis across
    // components would fail this and pass any tolerance written against pi.
    REQUIRE(recovered(1) == 0.0);
    REQUIRE(recovered(2) == 0.0);
}

TEST_CASE("SO3 very small rotation", "[so3][hardening][robustness]")
{
    ctrlpp::Vector<double, 3> phi;
    phi << 1e-12, 0.0, 0.0;

    auto q = ctrlpp::so3::exp(phi);
    REQUIRE_THAT(q.norm(), WithinAbs(1.0, so3_exp_norm_ops * so3_eps));

    // The scalar part is EXACTLY one: the cosine of half of 1e-12 differs from
    // one by about 1.25e-25, which is far below half an ulp at one, so the only
    // correctly rounded answer is one itself. A tolerance here would accept an
    // implementation whose small-angle handling was wrong by a hundred ulps.
    REQUIRE(q.w() == 1.0);

    auto recovered = ctrlpp::so3::log(q);

    // Exact, and exactly recoverable: the exponential's Taylor branch multiplies
    // the rotation vector by one half, which is a power of two and so exact, and
    // the logarithm's Taylor branch multiplies by two, which undoes it exactly.
    // An absolute tolerance at this scale is meaningless -- the value under test
    // is 1e-12 and the tolerance it replaces was a hundred times larger, so a
    // one percent error in the branch this case exists to check used to pass.
    REQUIRE(recovered(0) == 1e-12);
    REQUIRE(recovered(1) == 0.0);
    REQUIRE(recovered(2) == 0.0);
}

TEST_CASE("SO3 exp/log round trip across the small-angle branch threshold",
          "[so3][hardening][precision]")
{
    // The case above pins a value but NOT a branch: at an angle of 1e-12 the
    // Taylor branch and the arc-tangent branch agree to the last bit, so it
    // cannot tell which one ran. The exponential switches at an angle of 1e-7
    // and the logarithm at a vector norm of 1e-7, so the branch choice is only
    // observable within a factor of a few of that threshold. These angles
    // straddle it from both sides.
    const double budget = so3_exp_log_roundtrip_ops * so3_eps;

    for(const double theta : {1e-9, 1e-8, 9.9e-8, 1e-7, 1.1e-7, 1e-6, 1e-4, 1.0, 3.0})
    {
        CAPTURE(theta);
        ctrlpp::Vector<double, 3> phi;
        // An axis with two nonzero components, so a branch that mishandled the
        // vector norm could not hide behind a single-component rotation.
        phi << theta * 0.6, theta * 0.8, 0.0;

        auto const q = ctrlpp::so3::exp(phi);
        REQUIRE_THAT(q.norm(), WithinAbs(1.0, so3_exp_norm_ops * so3_eps));

        auto const recovered = ctrlpp::so3::log(q);
        // Relative to the angle, because the angles here span six decades and an
        // absolute budget would be vacuous at one end and unmeetable at the other.
        REQUIRE((recovered - phi).norm() <= budget * phi.norm());
    }
}
