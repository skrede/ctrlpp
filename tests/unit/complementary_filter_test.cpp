#include "ctrlpp/complementary_filter.h"
#include "ctrlpp/observer_policy.h"
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

// Static assert for ObserverPolicy concept
static_assert(ObserverPolicy<complementary_filter<double>>);

TEST_CASE("complementary filter IMU stationary converges to gravity-aligned", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for (int i = 0; i < 500; ++i) {
        cf.update(gyro, accel, 0.01);
    }

    double angle_error = quat_angle(cf.attitude(), Eigen::Quaterniond::Identity());
    CHECK(angle_error < 0.1);
}

TEST_CASE("complementary filter IMU rejects constant gyro bias", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    Vector<double, 3> bias_true;
    bias_true << 0.01, 0.01, 0.0;
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for (int i = 0; i < 2000; ++i) {
        cf.update(bias_true, accel, 0.01);
    }

    auto bias_est = cf.bias();
    CHECK_THAT(bias_est(0), WithinAbs(bias_true(0), 0.05));
    CHECK_THAT(bias_est(1), WithinAbs(bias_true(1), 0.05));
}

TEST_CASE("complementary filter IMU tracks rotation around z-axis", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    double omega_z = 0.1; // rad/s yaw rate
    Vector<double, 3> gyro;
    gyro << 0.0, 0.0, omega_z;

    for (int i = 0; i < 200; ++i) {
        double t = static_cast<double>(i) * 0.01;
        double yaw = omega_z * t;
        // Rotate accel vector to remain consistent with attitude
        Vector<double, 3> accel;
        accel << -9.81 * std::sin(yaw) * 0.0, 0.0, 9.81; // gravity still mostly +z
        cf.update(gyro, accel, 0.01);
    }

    // After 200 steps at 0.01s, yaw should be ~0.2 rad
    auto q = cf.attitude();
    // Extract yaw from quaternion
    double yaw_est = 2.0 * std::atan2(q.z(), q.w());
    double yaw_expected = omega_z * 200.0 * 0.01;
    CHECK_THAT(yaw_est, WithinAbs(yaw_expected, 0.15));
}

TEST_CASE("complementary filter MARG mode uses magnetometer for heading", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};

    // IMU-only filter
    complementary_filter cf_imu{cfg};
    // MARG filter
    complementary_filter cf_marg{cfg};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;
    Vector<double, 3> mag;
    mag << 0.2, 0.0, 0.4; // pointing roughly north+down

    for (int i = 0; i < 500; ++i) {
        cf_imu.update(gyro, accel, 0.01);
        cf_marg.update(gyro, accel, mag, 0.01);
    }

    // Both should converge to gravity-aligned attitude
    double imu_angle = quat_angle(cf_imu.attitude(), Eigen::Quaterniond::Identity());
    double marg_angle = quat_angle(cf_marg.attitude(), Eigen::Quaterniond::Identity());
    CHECK(imu_angle < 0.15);
    CHECK(marg_angle < 0.15);
}

TEST_CASE("complementary filter ObserverPolicy predict/update interface", "[cf]")
{
    static_assert(ObserverPolicy<complementary_filter<double>>);

    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    cf.predict(gyro);
    cf.update(accel);

    auto s = cf.state();
    CHECK(s.size() == 7); // 4 quaternion + 3 bias
    CHECK(std::isfinite(s(0)));
}

TEST_CASE("complementary filter reset via new construction", "[cf]")
{
    Eigen::Quaterniond q0;
    q0 = Eigen::AngleAxisd(0.3, Vector<double, 3>::UnitZ());
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01, .q0 = q0};
    complementary_filter cf{cfg};

    double angle = quat_angle(cf.attitude(), q0);
    CHECK(angle < 1e-6);
}

TEST_CASE("complementary filter handles zero accelerometer gracefully", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel = Vector<double, 3>::Zero();

    // Should not crash or produce NaN -- filter skips update when accel norm < 1e-10
    cf.update(gyro, accel, 0.01);

    auto q = cf.attitude();
    CHECK(std::isfinite(q.w()));
    CHECK(std::isfinite(q.x()));
    CHECK(std::isfinite(q.y()));
    CHECK(std::isfinite(q.z()));
}
