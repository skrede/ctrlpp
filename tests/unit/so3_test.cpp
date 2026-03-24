#include "ctrlpp/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <numbers>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace {

auto quat_angle(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) -> double
{
    return 2.0 * std::acos(std::min(1.0, std::abs(q1.dot(q2))));
}

} // namespace

TEST_CASE("so3 exp of zero vector returns identity quaternion", "[so3]")
{
    Vector<double, 3> zero = Vector<double, 3>::Zero();
    auto q = so3::exp(zero);
    CHECK_THAT(q.w(), WithinAbs(1.0, 1e-15));
    CHECK(q.vec().norm() < 1e-15);
}

TEST_CASE("so3 exp small angle uses Taylor approximation", "[so3]")
{
    Vector<double, 3> phi;
    phi << 1e-9, 0.0, 0.0;
    auto q = so3::exp(phi);
    CHECK_THAT(q.norm(), WithinAbs(1.0, 1e-14));
}

TEST_CASE("so3 exp at pi produces valid unit quaternion", "[so3]")
{
    Vector<double, 3> phi;
    phi << std::numbers::pi, 0.0, 0.0;
    auto q = so3::exp(phi);
    CHECK_THAT(q.norm(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("so3 log of identity returns zero vector", "[so3]")
{
    auto v = so3::log(Eigen::Quaterniond::Identity());
    CHECK(v.norm() < 1e-15);
}

TEST_CASE("so3 exp/log roundtrip for moderate angle", "[so3]")
{
    Vector<double, 3> phi;
    phi << 0.5, 0.3, -0.7;
    auto q = so3::exp(phi);
    auto recovered = so3::log(q);
    CHECK_THAT(recovered(0), WithinAbs(phi(0), 1e-12));
    CHECK_THAT(recovered(1), WithinAbs(phi(1), 1e-12));
    CHECK_THAT(recovered(2), WithinAbs(phi(2), 1e-12));
}

TEST_CASE("so3 compose is quaternion multiplication", "[so3]")
{
    Vector<double, 3> phi1, phi2;
    phi1 << 0.3, 0.1, -0.2;
    phi2 << -0.1, 0.4, 0.15;
    auto q1 = so3::exp(phi1);
    auto q2 = so3::exp(phi2);

    auto composed = so3::compose(q1, q2);
    Eigen::Quaterniond direct = q1 * q2;

    // Verify via rotation matrices
    auto R_composed = composed.toRotationMatrix();
    auto R_direct = direct.toRotationMatrix();
    CHECK((R_composed - R_direct).norm() < 1e-14);
}

TEST_CASE("so3 conjugate inverts rotation", "[so3]")
{
    Vector<double, 3> phi;
    phi << 0.5, -0.3, 0.7;
    auto q = so3::exp(phi);
    auto q_inv = so3::conjugate(q);
    auto identity = so3::compose(q, q_inv);

    CHECK_THAT(identity.w(), WithinAbs(1.0, 1e-14));
    CHECK(identity.vec().norm() < 1e-14);
}

TEST_CASE("so3 normalize produces unit quaternion", "[so3]")
{
    Eigen::Quaterniond q;
    q.w() = 2.0;
    q.x() = 1.0;
    q.y() = 0.0;
    q.z() = 0.0;
    auto qn = so3::normalize(q);
    CHECK_THAT(qn.norm(), WithinAbs(1.0, 1e-15));
}

TEST_CASE("so3 skew antisymmetry", "[so3]")
{
    Vector<double, 3> v;
    v << 1.5, -2.3, 0.7;
    auto S = so3::skew(v);
    auto diff = (S + S.transpose()).eval();
    CHECK(diff.norm() < 1e-15);
}

TEST_CASE("so3 to_vec/from_vec roundtrip", "[so3]")
{
    Vector<double, 3> phi;
    phi << 0.4, -0.2, 0.6;
    auto q = so3::exp(phi);
    auto v = so3::to_vec(q);
    auto q2 = so3::from_vec(v);

    CHECK_THAT(q2.w(), WithinAbs(q.w(), 1e-15));
    CHECK_THAT(q2.x(), WithinAbs(q.x(), 1e-15));
    CHECK_THAT(q2.y(), WithinAbs(q.y(), 1e-15));
    CHECK_THAT(q2.z(), WithinAbs(q.z(), 1e-15));
}

TEST_CASE("so3 log canonicalizes antipodal quaternion", "[so3]")
{
    Vector<double, 3> phi;
    phi << 0.5, 0.3, -0.7;
    auto q = so3::exp(phi);

    // Negate: antipodal quaternion represents same rotation
    Eigen::Quaterniond q_neg;
    q_neg.w() = -q.w();
    q_neg.x() = -q.x();
    q_neg.y() = -q.y();
    q_neg.z() = -q.z();

    auto log_pos = so3::log(q);
    auto log_neg = so3::log(q_neg);

    CHECK_THAT(log_neg(0), WithinAbs(log_pos(0), 1e-12));
    CHECK_THAT(log_neg(1), WithinAbs(log_pos(1), 1e-12));
    CHECK_THAT(log_neg(2), WithinAbs(log_pos(2), 1e-12));
}
