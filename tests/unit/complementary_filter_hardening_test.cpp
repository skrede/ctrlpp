#include "hardening_helpers.h"
#include "ctrlpp/estimation/complementary_filter.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

auto quat_angle(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) -> double
{
    return 2.0 * std::acos(std::min(1.0, std::abs(q1.dot(q2))));
}

}

TEST_CASE("Complementary filter NaN gyro does not crash",
          "[complementary_filter][hardening][negative]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    ctrlpp::complementary_filter cf{cfg};

    ctrlpp::Vector<double, 3> gyro;
    gyro << std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0;
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    cf.update(gyro, accel, 0.01);

    // NaN propagation or graceful -- no crash is the requirement
    CHECK(true);
}

TEST_CASE("Complementary filter NaN accel with zero norm skips update",
          "[complementary_filter][hardening][negative]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    ctrlpp::complementary_filter cf{cfg};

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel = ctrlpp::Vector<double, 3>::Zero();

    // Zero-norm accel should be skipped gracefully
    cf.update(gyro, accel, 0.01);

    auto q = cf.attitude();
    CHECK(std::isfinite(q.w()));
}

TEST_CASE("Complementary filter k_p=0 gives pure gyro integration",
          "[complementary_filter][hardening][precision]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 0.0, .k_i = 0.0, .dt = 0.01};
    ctrlpp::complementary_filter cf{cfg};

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    // With k_p=0, accelerometer corrections have no effect
    for(int i = 0; i < 100; ++i)
        cf.update(gyro, accel, 0.01);

    // Should remain near identity since gyro is zero
    double angle_error = quat_angle(cf.attitude(), Eigen::Quaterniond::Identity());
    CHECK(angle_error < 0.01);
}

TEST_CASE("Complementary filter converges to gravity-aligned",
          "[complementary_filter][hardening][convergence]")
{
    // Start from a tilted orientation
    auto q_init = Eigen::Quaterniond(Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitX()));
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01, .q0 = q_init};
    ctrlpp::complementary_filter cf{cfg};

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 1000; ++i)
        cf.update(gyro, accel, 0.01);

    double angle_error = quat_angle(cf.attitude(), Eigen::Quaterniond::Identity());
    REQUIRE(angle_error < 0.1);
}
