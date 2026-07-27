#include "ctrlpp/complementary_filter.h"
#include "ctrlpp/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace {

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. The rejection
// cases below do not use this helper: they assert the specific enumerator.
auto make_filter(const cf_config<double>& cfg) -> complementary_filter<double>
{
    auto created = complementary_filter<double>::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

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
    auto cf = make_filter(cfg);

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
    auto cf = make_filter(cfg);

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

TEST_CASE("complementary filter IMU bias estimate converges over a long horizon", "[cf]")
{
    // The corrected Mahony PI observer drives the gyro-bias estimate to the true
    // constant bias. Linearizing the single-axis loop (tilt error theta, bias
    // error b) gives theta_dot = -k_p theta - b, b_dot = k_i theta, i.e. the
    // characteristic polynomial s^2 + k_p s + k_i whose slow mode has time
    // constant tau = k_p / k_i. Integrating for T = 5 tau collapses the initial
    // bias error b0 by exp(-T/tau); the tolerance is that analytic decay envelope
    // with headroom for the fast-mode transient and the sin-vs-linear tilt term,
    // rather than a loose fixed band. The 20 s test above cannot see this ~400 s
    // time constant.
    const double k_p = 2.0, k_i = 0.005, dt = 0.01;
    cf_config<double> cfg{.k_p = k_p, .k_i = k_i, .dt = dt};
    auto cf = make_filter(cfg);

    Vector<double, 3> bias_true;
    bias_true << 0.01, 0.01, 0.0; // roll/pitch bias is observable from gravity; yaw is not
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    const double tau = k_p / k_i;      // slow-mode time constant of the PI loop
    const double horizon = 5.0 * tau;  // >= 500 s (here 2000 s)
    const int steps = static_cast<int>(horizon / dt);
    for (int i = 0; i < steps; ++i)
        cf.update(bias_true, accel, dt);

    const double b0 = bias_true(0); // initial bias error: the estimate starts at zero
    const double envelope = b0 * std::exp(-horizon / tau);
    // The measured residual tracks this slow-mode envelope closely (the fast mode
    // has fully decayed and the slow-mode participation is just under one), so a
    // 2x band asserts convergence to the analytic decay rather than a loose fixed
    // tolerance.
    const double tol = 2.0 * envelope;
    auto bias_est = cf.bias();
    INFO("bias(0) residual = " << std::abs(bias_est(0) - bias_true(0)) << ", envelope = " << envelope);
    INFO("bias(1) residual = " << std::abs(bias_est(1) - bias_true(1)) << ", envelope = " << envelope);
    CHECK_THAT(bias_est(0), WithinAbs(bias_true(0), tol));
    CHECK_THAT(bias_est(1), WithinAbs(bias_true(1), tol));
}

TEST_CASE("complementary filter IMU tracks rotation around z-axis", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    auto cf = make_filter(cfg);

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

    // Rotation is about the gravity axis, so the estimated gravity direction in
    // the body frame is invariant and the accel correction e vanishes exactly.
    // The filter therefore integrates the gyro at full rate about a fixed axis,
    // which composes exactly: yaw = omega_z * 200 * 0.01 = 0.2 rad, up to the
    // rounding accumulated across the quaternion normalizations.
    auto q = cf.attitude();
    // Extract yaw from quaternion
    double yaw_est = 2.0 * std::atan2(q.z(), q.w());
    double yaw_expected = omega_z * 200.0 * 0.01;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double integ_tol = 200.0 * eps * 8.0; // O(steps) normalizations, few-ulp headroom
    CHECK_THAT(yaw_est, WithinAbs(yaw_expected, integ_tol));
}

TEST_CASE("complementary filter MARG mode uses magnetometer for heading", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};

    // IMU-only filter
    auto cf_imu = make_filter(cfg);
    // MARG filter
    auto cf_marg = make_filter(cfg);

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
    auto cf = make_filter(cfg);

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
    auto cf = make_filter(cfg);

    double angle = quat_angle(cf.attitude(), q0);
    CHECK(angle < 1e-6);
}

TEST_CASE("complementary filter handles zero accelerometer gracefully", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    auto cf = make_filter(cfg);

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
    auto cf = make_filter(cfg);

    // With no proportional or integral gain, only gyro drives the filter
    Vector<double, 3> gyro;
    gyro << 0.0, 0.0, 0.1;
    Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 100; ++i)
        cf.update(gyro, accel, 0.01);

    // With k_p=0 and k_i=0 the filter is a pure gyro integrator. A fixed-axis
    // rotation composes exactly, so yaw = omega_z * steps * dt = 0.1 * 100 * 0.01
    // = 0.1 rad, up to the rounding accumulated across the quaternion
    // normalizations. The exponential map already carries the half-angle, so the
    // integrated angle is the full rate * dt (no extra 0.5 factor).
    auto q = cf.attitude();
    double yaw = 2.0 * std::atan2(q.z(), q.w());
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double yaw_expected = 0.1 * 100.0 * 0.01;
    const double integ_tol = 100.0 * eps * 8.0; // O(steps) normalizations, few-ulp headroom
    CHECK_THAT(yaw, WithinAbs(yaw_expected, integ_tol));
}

TEST_CASE("complementary filter k_i=0 has no bias estimation", "[cf]")
{
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.0, .dt = 0.01};
    auto cf = make_filter(cfg);

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
    auto cf_marg = make_filter(cfg);
    auto cf_imu = make_filter(cfg);

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
    auto cf = make_filter(cfg);

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
    auto cf = make_filter(cfg);

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
    auto cf = make_filter(cfg);

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
    auto cf = make_filter(cfg);

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
    // Start both from a tilted initial orientation
    Eigen::Quaterniond q0 = Eigen::Quaterniond(Eigen::AngleAxisd(0.5, Vector<double, 3>::UnitX()));
    cf_config<double> cfg_high_tilted{.k_p = 20.0, .k_i = 0.0, .dt = 0.01, .q0 = q0};
    cf_config<double> cfg_low_tilted{.k_p = 0.5, .k_i = 0.0, .dt = 0.01, .q0 = q0};
    auto cf_h = make_filter(cfg_high_tilted);
    auto cf_l = make_filter(cfg_low_tilted);

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

TEST_CASE("complementary filter create rejects a zero initial quaternion", "[cf]")
{
    // Before the fallible factory existed, the constructor stored the zero
    // quaternion raw; the correction terms then used it as a rotation and the
    // attitude estimate degenerated silently.
    Eigen::Quaterniond q0{0.0, 0.0, 0.0, 0.0};
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01, .q0 = q0};

    auto cf = complementary_filter<double>::create(cfg);
    REQUIRE_FALSE(cf.has_value());
    CHECK(cf.error() == filter_error::degenerate_quaternion);
}

TEST_CASE("complementary filter create accepts and normalizes a non-unit quaternion", "[cf]")
{
    // The guard rejects only degeneracy: any finite nonzero quaternion is
    // accepted and brought onto the unit sphere at construction, since the
    // correction terms treat the stored quaternion as a unit rotation.
    Eigen::Quaterniond q0{2.0, 0.0, 0.0, 0.0};
    cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01, .q0 = q0};

    auto cf = complementary_filter<double>::create(cfg);
    REQUIRE(cf.has_value());
    CHECK((cf->attitude().coeffs().array() == q0.normalized().coeffs().array()).all());
}
