#include "hardening_helpers.h"
#include "ctrlpp/estimation/complementary_filter.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace {

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here.
auto make_filter(const ctrlpp::cf_config<double>& cfg) -> ctrlpp::complementary_filter<double>
{
    auto created = ctrlpp::complementary_filter<double>::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

auto quat_angle(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) -> double
{
    return 2.0 * std::acos(std::min(1.0, std::abs(q1.dot(q2))));
}

}

TEST_CASE("Complementary filter NaN gyro is rejected without touching the attitude",
          "[complementary_filter][hardening][negative]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    auto cf = make_filter(cfg);
    // Stepped only with the valid rate, never with the poisoned one, so it says
    // what the filter would have carried had the bad sample never arrived.
    auto reference = make_filter(cfg);

    ctrlpp::Vector<double, 3> gyro_bad;
    gyro_bad << std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0;
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    // Snapshot immediately before the poisoned step. The filter carries no
    // covariance, so the preserved-invariant argument has an attitude and bias
    // half only.
    const Eigen::Quaterniond q_before = cf.attitude();
    const ctrlpp::Vector<double, 3> bias_before = cf.bias();

    const auto rejected = cf.update(gyro_bad, accel, 0.01);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::cf_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried attitude at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(cf.attitude().coeffs() == q_before.coeffs());
    CHECK(cf.bias() == bias_before);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(cf.health() == ctrlpp::cf_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    const ctrlpp::Vector<double, 3> gyro_good = ctrlpp::Vector<double, 3>::Zero();
    REQUIRE(cf.update(gyro_good, accel, 0.01).has_value());
    REQUIRE(reference.update(gyro_good, accel, 0.01).has_value());

    CHECK(cf.attitude().coeffs() == reference.attitude().coeffs());
    CHECK(cf.bias() == reference.bias());
}

TEST_CASE("Complementary filter non-finite timestep is rejected and names the clock",
          "[complementary_filter][hardening][negative]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    auto cf = make_filter(cfg);

    const ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    const Eigen::Quaterniond q_before = cf.attitude();

    // A broken clock and a broken sensor are different faults with different
    // repairs, so they carry different enumerators.
    const auto rejected = cf.update(gyro, accel, std::numeric_limits<double>::quiet_NaN());

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::cf_update_error::non_finite_timestep);
    CHECK(cf.attitude().coeffs() == q_before.coeffs());
}

TEST_CASE("Complementary filter zero-norm accel skips the correction and still succeeds",
          "[complementary_filter][hardening][negative]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};
    auto cf = make_filter(cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel = ctrlpp::Vector<double, 3>::Zero();

    // A zero-norm acceleration carries no gravity direction, so there is no
    // correction to apply. That is a SUCCESS in which the commanded correction
    // was not performed, not a failure: the step did exactly what the algorithm
    // prescribes. Reporting it on the failure channel would teach the caller
    // that channel carries non-failures.
    REQUIRE(cf.update(gyro, accel, 0.01).has_value());

    auto q = cf.attitude();
    CHECK(std::isfinite(q.w()));
}

TEST_CASE("Complementary filter k_p=0 gives pure gyro integration",
          "[complementary_filter][hardening][precision]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 0.0, .k_i = 0.0, .dt = 0.01};
    auto cf = make_filter(cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    // With k_p=0, accelerometer corrections have no effect
    for(int i = 0; i < 100; ++i)
        REQUIRE(cf.update(gyro, accel, 0.01).has_value());

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
    auto cf = make_filter(cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 1000; ++i)
        REQUIRE(cf.update(gyro, accel, 0.01).has_value());

    double angle_error = quat_angle(cf.attitude(), Eigen::Quaterniond::Identity());
    REQUIRE(angle_error < 0.1);
}
