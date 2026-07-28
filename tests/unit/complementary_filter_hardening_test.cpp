// What the oracles in this file decide.
//
// The filter carries an attitude quaternion and a gyro bias and nothing else,
// so every claim here is a claim about those two:
//
//  * A rejected step leaves both BITWISE unchanged, names the specific cause,
//    and does not degrade the health status, and a later valid step produces
//    exactly what an instance that never saw the bad sample produces.
//  * A step whose accelerometer carries no gravity direction skips the
//    correction. With a zero rate there is then no arithmetic left to perform,
//    so the attitude and the bias are BITWISE what they were. That is an exact
//    statement and a tolerance would admit a step that partially ran.
//  * With both correction gains and the rate at zero, every increment is an
//    exact product of zero and the quaternion composition is with the exact
//    identity, so after a hundred steps the attitude is BITWISE the identity.
//    Measured: it is, and the per-step renormalization introduces nothing,
//    because the norm it divides by is exactly one.
//  * Convergence is asserted against the CLOSED-LOOP MAP the configuration
//    fixes, not against a decay envelope. With the rate at zero and the
//    measured gravity exactly that of the identity attitude, the whole run
//    stays in the one-parameter subgroup about the tilt axis, where the filter
//    reduces exactly to
//        e     = -sin(theta)
//        b    <- b - k_i e dt
//        theta<- theta + dt (-b + k_p e)
//    with e the gravity cross product, b the estimated bias and theta the
//    signed tilt. That recursion is carried alongside the filter from the same
//    configured gains and the realized tilt is asserted against it.
//
// What they deliberately do not decide. Nothing here asserts a decay RATE. The
// proportional-integral pair places two closed-loop roots whose slow one is
// three orders below the fast one, so an envelope written from the proportional
// gain alone is wrong by six orders at this horizon -- measured: the realized
// residual after ten seconds is 6.1e-4 radians, not the 1e-9 such an envelope
// predicts, because the integral term parks a bias of 1.2e-3 and the tilt
// settles at that bias over the proportional gain. The map above states the
// whole behavior and needs no rate.

#include "hardening_helpers.h"
#include "ctrlpp/estimation/complementary_filter.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

namespace {

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here.
auto make_filter(const ctrlpp::cf_config<double>& cfg) -> ctrlpp::complementary_filter<double>
{
    auto created = ctrlpp::complementary_filter<double>::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

// Signed tilt about the x axis, read straight off the quaternion rather than
// through an arc cosine. The arc cosine of a dot product loses half its digits
// near a zero angle -- it saturates to exactly zero once the vector part falls
// below about 1e-8 -- and the convergence case below tracks a tilt four orders
// under that.
auto tilt_x(const Eigen::Quaterniond& q) -> double
{
    return 2.0 * std::atan2(q.x(), q.w());
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
    const Eigen::Quaterniond q_before = cf.attitude();
    const ctrlpp::Vector<double, 3> bias_before = cf.bias();

    REQUIRE(cf.update(gyro, accel, 0.01).has_value());

    // With no correction AND a zero rate there is nothing left for the step to
    // compute, so the carried estimate is BITWISE what it was. Exact, because
    // the skip performs no arithmetic on it -- the same argument the rejection
    // cases above make, on a path that succeeds.
    CHECK(cf.attitude().coeffs() == q_before.coeffs());
    CHECK(cf.bias() == bias_before);
    // Skipping a correction is not damage: the estimate is exactly as good as
    // it was before the sample arrived, so the status must still say so.
    CHECK(cf.health() == ctrlpp::cf_health::ok);
}

TEST_CASE("Complementary filter k_p=0 gives pure gyro integration",
          "[complementary_filter][hardening][precision]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 0.0, .k_i = 0.0, .dt = 0.01};
    auto cf = make_filter(cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    for(int i = 0; i < 100; ++i)
        REQUIRE(cf.update(gyro, accel, 0.01).has_value());

    // Every term of the update is an exact product of zero: both correction
    // gains are zero, the rate is zero, so the tangent vector is exactly the
    // zero vector, the exponential of it is exactly the identity quaternion,
    // the composition with the identity reproduces every coefficient, and the
    // renormalization divides by a norm of exactly one. A hundred repetitions
    // of an exact identity is still an identity, so this is bitwise, not near.
    CHECK(cf.attitude().coeffs() == Eigen::Quaterniond::Identity().coeffs());
    CHECK(cf.bias() == ctrlpp::Vector<double, 3>::Zero());
}

TEST_CASE("Complementary filter converges to gravity-aligned",
          "[complementary_filter][hardening][convergence]")
{
    // Start from a tilted orientation
    constexpr double k_p = 2.0;
    constexpr double k_i = 0.005;
    constexpr double step = 0.01;
    constexpr double initial_tilt = 0.5;
    constexpr int cycles = 1000;

    auto q_init = Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitX()));
    ctrlpp::cf_config<double> cfg{.k_p = k_p, .k_i = k_i, .dt = step, .q0 = q_init};
    auto cf = make_filter(cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    ctrlpp::Vector<double, 3> accel;
    accel << 0.0, 0.0, 9.81;

    // The closed-loop map, carried from the configured gains alone. For a tilt
    // theta about x the normalized acceleration is (0,0,1) and the estimated
    // gravity in the body frame is (0, sin theta, cos theta), so their cross
    // product is (-sin theta, 0, 0): the correction is exactly -sin(theta)
    // about the same axis. The bias moves first, then enters the corrected
    // rate, which is the order `integrate_gyro` uses.
    //
    // Thirty-five roundings enter one cycle: the sine (1), the bias product and
    // its subtraction (2, 3), the proportional product and the two sums forming
    // the corrected rate (4, 5, 6), the tangent scaling (7), the exponential's
    // half angle, cosine, sine and cardinal-sine scaling (8-11), the sixteen
    // products and sums of the quaternion composition (12-27), and the
    // renormalization's four squares, three sums, square root and division
    // (28-35). They do NOT accumulate over the run: the map is a contraction in
    // the tilt, so a rounding injected at one cycle is attenuated by every
    // cycle after it and the realized disagreement stays at the level of a
    // single cycle. The budget is therefore one cycle's roundings at the
    // initial tilt's scale, and the realized worst is 0.75 of one such
    // rounding.
    constexpr int cycle_rounding_ops = 35;
    const double budget = cycle_rounding_ops * std::numeric_limits<double>::epsilon() * initial_tilt;

    double theta = initial_tilt;
    double bias = 0.0;
    double worst = 0.0;

    for(int i = 0; i < cycles; ++i)
    {
        REQUIRE(cf.update(gyro, accel, step).has_value());

        const double correction = -std::sin(theta);
        bias -= k_i * correction * step;
        theta += step * (-bias + k_p * correction);

        worst = std::max(worst, std::abs(tilt_x(cf.attitude()) - theta));
        REQUIRE(std::abs(cf.bias()(0) - bias) <= budget);
    }

    CHECK(worst <= budget);

    // Convergence itself, stated without a threshold. The tilt ends smaller in
    // magnitude than it began and on the OTHER side of zero: only an integral
    // term can carry it past the setpoint, so a controller that lost its
    // integral action, or applied it with the wrong sign, fails this without
    // any number having to be chosen. The estimated bias ends positive for the
    // same reason -- it is the integral of a correction that was negative
    // throughout the approach.
    CHECK(std::abs(tilt_x(cf.attitude())) < initial_tilt);
    CHECK(tilt_x(cf.attitude()) < 0.0);
    CHECK(cf.bias()(0) > 0.0);
}
