#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>

#include <Eigen/Geometry>

#include <cmath>

using namespace ctrlpp;

namespace
{

auto gen_rotation_vector(double max_angle) -> rc::Gen<Vector<double, 3>>
{
    return rc::gen::map(rc::gen::tuple(rc::gen::inRange(-1000000, 1000000), rc::gen::inRange(-1000000, 1000000), rc::gen::inRange(-1000000, 1000000)),
                        [max_angle](const std::tuple<int, int, int>& t)
                        {
                            Vector<double, 3> v;
                            v(0) = max_angle * static_cast<double>(std::get<0>(t)) / 1000000.0;
                            v(1) = max_angle * static_cast<double>(std::get<1>(t)) / 1000000.0;
                            v(2) = max_angle * static_cast<double>(std::get<2>(t)) / 1000000.0;
                            // Scale to keep within max_angle
                            double n = v.norm();
                            if(n > max_angle)
                                v *= max_angle / n;
                            return v;
                        });
}

auto gen_unit_quaternion() -> rc::Gen<Eigen::Quaternion<double>>
{
    return rc::gen::map(gen_rotation_vector(3.0), [](const Vector<double, 3>& phi) { return so3::exp(phi); });
}

} // namespace

TEST_CASE("so3 property tests", "[so3][property]")
{
    SECTION("so3 exp/log roundtrip preserves rotation vector")
    {
        rc::prop("exp/log roundtrip", [](void)
                 {
            auto phi = *gen_rotation_vector(3.0);
            auto q = so3::exp(phi);
            auto phi_back = so3::log(q);
            auto q_back = so3::exp(phi_back);

            // quaternions represent same rotation: |q . q_back| > 1 - eps
            RC_ASSERT(std::abs(q.dot(q_back)) > 1.0 - 1e-10); });
    }

    SECTION("so3 compose is associative")
    {
        rc::prop("associativity", [](void)
                 {
            auto q1 = *gen_unit_quaternion();
            auto q2 = *gen_unit_quaternion();
            auto q3 = *gen_unit_quaternion();

            auto lhs = so3::compose(so3::compose(q1, q2), q3);
            auto rhs = so3::compose(q1, so3::compose(q2, q3));

            RC_ASSERT(std::abs(lhs.dot(rhs)) > 1.0 - 1e-10); });
    }

    SECTION("so3 identity is neutral element")
    {
        rc::prop("identity", [](void)
                 {
            auto q = *gen_unit_quaternion();
            Eigen::Quaternion<double> id = Eigen::Quaternion<double>::Identity();

            auto lq = so3::compose(id, q);
            auto rq = so3::compose(q, id);

            RC_ASSERT(std::abs(lq.dot(q)) > 1.0 - 1e-10);
            RC_ASSERT(std::abs(rq.dot(q)) > 1.0 - 1e-10); });
    }

    SECTION("so3 conjugate is group inverse")
    {
        rc::prop("inverse", [](void)
                 {
            auto q = *gen_unit_quaternion();
            auto qinv = so3::conjugate(q);
            auto result = so3::compose(q, qinv);

            Eigen::Quaternion<double> id = Eigen::Quaternion<double>::Identity();
            RC_ASSERT(std::abs(result.dot(id)) > 1.0 - 1e-10); });
    }

    SECTION("so3 normalize is idempotent")
    {
        rc::prop("normalize idempotent", [](void)
                 {
            auto q = *gen_unit_quaternion();
            auto qn = so3::normalize(q);
            auto qnn = so3::normalize(qn);

            RC_ASSERT(std::abs(qn.w() - qnn.w()) < 1e-15);
            RC_ASSERT((qn.vec() - qnn.vec()).norm() < 1e-15); });
    }

    SECTION("so3 left Jacobian is consistent with finite-difference approximation")
    {
        rc::prop("Jacobian consistency", [](void)
                 {
            auto phi = *gen_rotation_vector(2.5); // stay away from pi
            double theta = phi.norm();

            // Compute analytical left Jacobian J_l(phi)
            // Sola2018 Eq. 145: J_l = sin(t)/t * I + (1 - sin(t)/t) * a*a^T + (1 - cos(t))/t * [a]_x
            Eigen::Matrix3d J_l;
            if(theta < 1e-7)
            {
                // Taylor: J_l -> I + 0.5 * [phi]_x
                J_l = Eigen::Matrix3d::Identity() + 0.5 * so3::skew(phi);
            }
            else
            {
                Vector<double, 3> a = phi / theta;
                double st = std::sin(theta);
                double ct = std::cos(theta);
                J_l = (st / theta) * Eigen::Matrix3d::Identity() + (1.0 - st / theta) * a * a.transpose() + ((1.0 - ct) / theta) * so3::skew(a);
            }

            // Numerical left Jacobian: column j = log(exp(phi + eps*e_j) * exp(phi)^{-1}) / eps
            constexpr double eps = 1e-7;
            Eigen::Matrix3d J_num;
            auto q_base = so3::exp(phi);

            for(int j = 0; j < 3; ++j)
            {
                Vector<double, 3> e_j = Vector<double, 3>::Zero();
                e_j(j) = eps;
                Vector<double, 3> phi_p = (phi + e_j).eval();
                auto q_p = so3::exp(phi_p);
                auto delta = so3::compose(q_p, so3::conjugate(q_base));
                J_num.col(j) = so3::log(delta) / eps;
            }

            double err = (J_l - J_num).norm();
            RC_ASSERT(err < 1e-4); });
    }
}
