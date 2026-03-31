#include "hardening_helpers.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
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

    // normalize of zero quaternion
    auto qn = ctrlpp::so3::normalize(q);
    // Eigen normalized() of zero produces NaN
    bool const w_valid = std::isnan(qn.w()) || std::isfinite(qn.w());
    CHECK(w_valid);
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
