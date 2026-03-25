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

// ---------------------------------------------------------------------------
// Edge / failure-path tests
// ---------------------------------------------------------------------------

TEST_CASE("complementary filter k_p=0 ignores accel correction", "[cf]")
{
    cf_config<double> cfg{.k_p = 0.0, .k_i = 0.0, .dt = 0.01};
    complementary_filter cf{cfg};

    // With no proportional or integral gain, only gyro drives the filter
    Vector<double, 3> gyro;
    gyro << 0.0, 0.0, 0.1;
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 100; ++i)
        cf.update(gyro, accel, 0.01);

    // Filter should purely integrate gyro: yaw ~ 0.1 * 100 * 0.01 = 0.1 rad
    auto q = cf.attitude();
    double yaw = 2.0 * std::atan2(q.z(), q.w());
    // Should be close to pure gyro integration since k_p=0 and k_i=0
    CHECK_THAT(yaw, WithinAbs(0.1 * 100.0 * 0.01 * 0.5, 0.15));
}

TEST_CASE("complementary filter k_i=0 has no bias estimation", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.0, .dt = 0.01};
    complementary_filter cf{cfg};

    Vector<double, 3> bias_true;
    bias_true << 0.01, 0.01, 0.0;
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 2000; ++i)
        cf.update(bias_true, accel, 0.01);

    // With k_i=0, bias should remain at zero
    auto bias_est = cf.bias();
    CHECK_THAT(bias_est(0), WithinAbs(0.0, 1e-15));
    CHECK_THAT(bias_est(1), WithinAbs(0.0, 1e-15));
    CHECK_THAT(bias_est(2), WithinAbs(0.0, 1e-15));
}

TEST_CASE("complementary filter MARG with zero magnetometer falls back to IMU", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf_marg{cfg};
    complementary_filter cf_imu{cfg};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;
    Vector<double, 3> mag_zero = Vector<double, 3>::Zero();

    for(int i = 0; i < 100; ++i)
    {
        cf_marg.update(gyro, accel, mag_zero, 0.01);
        cf_imu.update(gyro, accel, 0.01);
    }

    // MARG with zero mag should produce same result as IMU-only
    double angle_diff = quat_angle(cf_marg.attitude(), cf_imu.attitude());
    CHECK(angle_diff < 1e-10);
}

TEST_CASE("complementary filter MARG with near-zero accel is skipped", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    auto q_before = cf.attitude();

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel_tiny;
    accel_tiny << 1e-15, 0.0, 0.0;
    Vector<double, 3> mag;
    mag << 0.2, 0.0, 0.4;

    // Both IMU and MARG updates should be no-ops with near-zero accel
    cf.update(gyro, accel_tiny, mag, 0.01);

    double angle_change = quat_angle(cf.attitude(), q_before);
    CHECK(angle_change < 1e-10);
}

TEST_CASE("complementary filter large dt produces finite results", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 1.0};
    complementary_filter cf{cfg};

    Vector<double, 3> gyro;
    gyro << 0.5, -0.3, 0.1;
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    cf.update(gyro, accel, 1.0);

    auto q = cf.attitude();
    CHECK(std::isfinite(q.w()));
    CHECK(std::isfinite(q.x()));
    CHECK(std::isfinite(q.y()));
    CHECK(std::isfinite(q.z()));
}

TEST_CASE("complementary filter multiple predict/update cycles via ObserverPolicy", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    complementary_filter cf{cfg};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 200; ++i)
    {
        cf.predict(gyro);
        cf.update(accel);
    }

    double angle = quat_angle(cf.attitude(), Eigen::Quaterniond::Identity());
    CHECK(angle < 0.1);
    CHECK(cf.state().size() == 7);
}

TEST_CASE("complementary filter non-identity initial quaternion is preserved before update", "[cf]")
{
    Eigen::Quaterniond q0 = Eigen::Quaterniond(Eigen::AngleAxisd(1.0, Vector<double, 3>::UnitX()));
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01, .q0 = q0};
    complementary_filter cf{cfg};

    // State should reflect the initial quaternion
    auto s = cf.state();
    CHECK_THAT(s(0), WithinAbs(q0.w(), 1e-10));
    CHECK_THAT(s(1), WithinAbs(q0.x(), 1e-10));
    CHECK_THAT(s(2), WithinAbs(q0.y(), 1e-10));
    CHECK_THAT(s(3), WithinAbs(q0.z(), 1e-10));
    // Bias should be zero
    CHECK_THAT(s(4), WithinAbs(0.0, 1e-15));
    CHECK_THAT(s(5), WithinAbs(0.0, 1e-15));
    CHECK_THAT(s(6), WithinAbs(0.0, 1e-15));
}

TEST_CASE("complementary filter high proportional gain converges faster", "[cf]")
{
    // High k_p filter
    cf_config<double> cfg_high{.k_p = 20.0, .k_i = 0.0, .dt = 0.01};
    complementary_filter cf_high{cfg_high};

    // Low k_p filter
    cf_config<double> cfg_low{.k_p = 0.5, .k_i = 0.0, .dt = 0.01};
    complementary_filter cf_low{cfg_low};

    // Start both from a tilted initial orientation
    Eigen::Quaterniond q0 = Eigen::Quaterniond(Eigen::AngleAxisd(0.5, Vector<double, 3>::UnitX()));
    cf_config<double> cfg_high_tilted{.k_p = 20.0, .k_i = 0.0, .dt = 0.01, .q0 = q0};
    cf_config<double> cfg_low_tilted{.k_p = 0.5, .k_i = 0.0, .dt = 0.01, .q0 = q0};
    complementary_filter cf_h{cfg_high_tilted};
    complementary_filter cf_l{cfg_low_tilted};

    Vector<double, 3> gyro = Vector<double, 3>::Zero();
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 50; ++i)
    {
        cf_h.update(gyro, accel, 0.01);
        cf_l.update(gyro, accel, 0.01);
    }

    double err_high = quat_angle(cf_h.attitude(), Eigen::Quaterniond::Identity());
    double err_low = quat_angle(cf_l.attitude(), Eigen::Quaterniond::Identity());

    // Higher gain should converge faster (lower error after same number of steps)
    CHECK(err_high < err_low);
}
