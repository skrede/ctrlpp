#include "hardening_helpers.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>

using Catch::Matchers::WithinAbs;

// ── SO(3) hardening ────────────────────────────────────────────────────────────

TEST_CASE("SO3 zero quaternion", "[so3][hardening][negative]")
{
    Eigen::Quaterniond q;
    q.w() = 0.0;
    q.x() = 0.0;
    q.y() = 0.0;
    q.z() = 0.0;

    // log of zero quaternion -- degenerate, should not crash
    auto v = ctrlpp::so3::log(q);
    // Result is implementation-defined but must not crash
    (void)v;

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

    // log should still work (just not unit)
    auto v = ctrlpp::so3::log(q);
    REQUIRE(std::isfinite(v(0)));
    REQUIRE(std::isfinite(v(1)));
    REQUIRE(std::isfinite(v(2)));
}

TEST_CASE("SO3 90-degree rotation around z-axis", "[so3][hardening][precision]")
{
    ctrlpp::Vector<double, 3> phi;
    phi << 0.0, 0.0, std::numbers::pi / 2.0;

    auto q = ctrlpp::so3::exp(phi);
    auto R = q.toRotationMatrix();

    // Expected: [[0,-1,0],[1,0,0],[0,0,1]]
    REQUIRE_THAT(R(0, 0), WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(R(0, 1), WithinAbs(-1.0, 1e-10));
    REQUIRE_THAT(R(1, 0), WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(R(1, 1), WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(R(2, 2), WithinAbs(1.0, 1e-10));
}

TEST_CASE("SO3 180-degree rotation (near singularity)", "[so3][hardening][robustness]")
{
    ctrlpp::Vector<double, 3> phi;
    phi << std::numbers::pi, 0.0, 0.0;

    auto q = ctrlpp::so3::exp(phi);
    REQUIRE_THAT(q.norm(), WithinAbs(1.0, 1e-10));

    // exp/log roundtrip -- 180 degrees is near the branch cut
    auto recovered = ctrlpp::so3::log(q);
    REQUIRE_THAT(recovered.norm(), WithinAbs(std::numbers::pi, 1e-8));
}

TEST_CASE("SO3 very small rotation", "[so3][hardening][robustness]")
{
    ctrlpp::Vector<double, 3> phi;
    phi << 1e-12, 0.0, 0.0;

    auto q = ctrlpp::so3::exp(phi);
    REQUIRE_THAT(q.norm(), WithinAbs(1.0, 1e-14));
    REQUIRE_THAT(q.w(), WithinAbs(1.0, 1e-14));

    auto recovered = ctrlpp::so3::log(q);
    REQUIRE_THAT(recovered(0), WithinAbs(1e-12, 1e-14));
    REQUIRE_THAT(recovered(1), WithinAbs(0.0, 1e-14));
    REQUIRE_THAT(recovered(2), WithinAbs(0.0, 1e-14));
}
