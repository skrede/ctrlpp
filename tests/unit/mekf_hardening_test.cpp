// What the oracles in this file decide.
//
// The filter carries a nominal attitude and bias plus a six-dimensional
// error-state covariance, and every case here is fed a gravity measurement.
// That measurement is BLIND TO ROTATION ABOUT GRAVITY, and this file's oracles
// are built around that fact rather than around finiteness:
//
//  * A rejected measurement leaves the estimate, the covariance and the
//    quaternion norm BITWISE unchanged, names its enumerator, and does not
//    degrade the health status.
//  * With no process noise, a zero rate and a measurement that is exactly the
//    gravity of the carried attitude, the innovation is EXACTLY zero, so the
//    correction is exactly zero and the nominal state is bitwise unchanged --
//    while the covariance still contracts in the two directions gravity
//    observes and does not move in the one it does not.
//  * Positive definiteness is asserted with NO negative floor. The symmetrized
//    recursion produces a covariance whose two triangles hold the same bits, so
//    symmetry is asserted EXACTLY, and the smallest eigenvalue is asserted
//    strictly positive rather than merely above a chosen negative number.
//  * Convergence is asserted three ways, none of them a fitted fraction of the
//    initial error: the true attitude is an EXACT FIXED POINT (an instance
//    seeded there never moves, bitwise, over the whole run); the tilt keeps
//    falling across four checkpoints an order apart, which no filter that
//    stalls at a floor can do; and the trajectory is BITWISE identical whether
//    the initial tilt is about x or about y, which is the rotational symmetry
//    of a gravity measurement and which an axis-asymmetric Jacobian breaks.
//  * The near-half-turn case is a case about UNOBSERVABILITY, not about a
//    branch cut. Its fixture puts the whole initial error in yaw, which
//    gravity cannot see, so the innovation is exactly zero at every step, the
//    estimate is BITWISE frozen, the two observed variances fall far below
//    their initial value and the unobserved one GROWS at every step. Asserting
//    convergence there would be red against entirely correct behavior.
//
// What they deliberately do not decide. No case asserts a convergence RATE. The
// rate is set by the filter's own Riccati recursion, and the only way to state
// it independently is to solve that recursion again in the test, which would
// rest the suite on a second implementation of the estimator under test. The
// fixed point, the no-floor ladder and the axis symmetry pin the behavior
// without it. The reported covariance is likewise NOT used as a convergence
// bound: measured, the realized error is four orders inside the filter's own
// one-sigma, so such a bound would be no tighter than the fitted number it
// would replace.

#include "hardening_helpers.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace {

struct gravity_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q,
                    const ctrlpp::Vector<double, 3>& /*b*/) const -> ctrlpp::Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

// The tilt angle, read off the quaternion's vector part rather than through an
// arc cosine of a dot product. The arc cosine saturates to exactly zero once
// the vector part falls below about 1e-8, and the convergence case tracks a
// tilt that passes through 1e-13.
auto rotation_angle(const Eigen::Quaterniond& q) -> double
{
    return 2.0 * std::atan2(q.vec().norm(), std::abs(q.w()));
}

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here.
auto build_mekf(const ctrlpp::mekf_config<double, 3, 3>& cfg)
    -> ctrlpp::mekf<double, 3, 3, gravity_measurement>
{
    auto created = ctrlpp::mekf<double, 3, 3, gravity_measurement>::create(gravity_measurement{}, cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

auto make_mekf()
{
    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.01;
    cfg.dt = 0.01;
    return build_mekf(cfg);
}

auto gravity_of_identity() -> ctrlpp::Vector<double, 3>
{
    ctrlpp::Vector<double, 3> z;
    z << 0.0, 0.0, 1.0;
    return z;
}

// A step that moves the attitude ends in a renormalization: one division of
// four coefficients by a norm, whose result is correctly rounded, so the
// realized norm is within one rounding of unity. The deviation does not
// accumulate across steps, because each step renormalizes from scratch.
// Measured worst over two thousand steps: half of one such rounding.
constexpr int normalize_rounding_ops = 1;
const double unit_norm_budget = normalize_rounding_ops * std::numeric_limits<double>::epsilon();

}

TEST_CASE("MEKF zero rotation noise covariance", "[mekf][hardening][negative]")
{
    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q = ctrlpp::Matrix<double, 6, 6>::Zero();
    cfg.R *= 0.01;
    cfg.dt = 0.01;

    auto filter = build_mekf(cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(gyro);

    const ctrlpp::Vector<double, 7> state_after_predict = filter.state();
    const ctrlpp::Matrix<double, 6, 6> P_after_predict = filter.covariance();

    REQUIRE(filter.update(gravity_of_identity()).has_value());

    // The measurement IS the gravity the carried attitude predicts, so the
    // innovation is a difference of identical vectors: exactly zero, not nearly
    // zero. The correction is the gain times that, so it is exactly the zero
    // vector and the nominal state comes through bitwise.
    CHECK(filter.innovation() == ctrlpp::Vector<double, 3>::Zero());
    CHECK(filter.state() == state_after_predict);
    CHECK(filter.attitude().norm() == 1.0);
    CHECK(filter.health() == ctrlpp::mekf_health::ok);

    // The covariance still moves, and it moves ANISOTROPICALLY, which is the
    // real content of a gravity update: the two attitude directions gravity
    // observes contract by two orders, while rotation ABOUT gravity is
    // unobserved and its variance comes through the update bitwise. With no
    // process noise there is nothing to reinflate it.
    CHECK(filter.covariance()(0, 0) < P_after_predict(0, 0));
    CHECK(filter.covariance()(1, 1) < P_after_predict(1, 1));
    CHECK(filter.covariance()(2, 2) == P_after_predict(2, 2));
    // The recursion symmetrizes, so the two triangles hold the same bits.
    CHECK(filter.covariance() == filter.covariance().transpose());
}

TEST_CASE("MEKF NaN quaternion measurement is rejected without touching the estimate",
          "[mekf][hardening][negative]")
{
    auto filter = make_mekf();
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = make_mekf();

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(gyro);
    reference.predict(gyro);

    // Snapshot immediately before the poisoned step.
    const ctrlpp::Vector<double, 7> state_before = filter.state();
    const ctrlpp::Matrix<double, 6, 6> P_before = filter.covariance();

    ctrlpp::Vector<double, 3> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0;

    const auto rejected = filter.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::mekf_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(filter.state() == state_before);
    // The covariance was already measurement-independent before the guard
    // existed -- update_covariance(K, H, delta_xi) takes the gain, the
    // measurement Jacobian and the correction, never z -- so this half of the
    // invariant is structural. The state half is what the guard adds.
    CHECK(filter.covariance() == P_before);
    // On a manifold the attitude quaternion's unit norm is an extra invariant,
    // and the rejection preserves it exactly rather than approximately.
    CHECK(filter.attitude().norm() == 1.0);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(filter.health() == ctrlpp::mekf_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    ctrlpp::Vector<double, 3> z_good;
    z_good << 0.0, 0.0, 1.0;
    REQUIRE(filter.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(filter.state() == reference.state());
    CHECK(filter.covariance() == reference.covariance());
}

TEST_CASE("MEKF covariance stays symmetric positive definite over 1000 steps",
          "[mekf][hardening][stability]")
{
    auto filter = make_mekf();

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();

    for(int k = 0; k < 1000; ++k)
    {
        filter.predict(gyro);
        REQUIRE(filter.update(gravity_of_identity()).has_value());

        // Symmetry is EXACT, not approximate: both the propagation and the
        // update pass their result through the symmetrizing helper, so the two
        // triangles hold the same bits. A tolerance here would accept a
        // recursion that had quietly stopped symmetrizing.
        REQUIRE(filter.covariance() == filter.covariance().transpose());

        // Positive definiteness with NO negative floor. A floor admits an
        // indefinite covariance in a case named for definiteness; the realized
        // smallest eigenvalue never approaches zero, so nothing has to be
        // allowed for.
        Eigen::SelfAdjointEigenSolver<ctrlpp::Matrix<double, 6, 6>> eigsolver(filter.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff() > 0.0);
    }
}

TEST_CASE("MEKF attitude converges to true orientation", "[mekf][hardening][convergence]")
{
    constexpr double initial_tilt = 0.3;

    auto configured = [](const Eigen::Quaterniond& q0) {
        ctrlpp::mekf_config<double, 3, 3> cfg;
        cfg.Q *= 1e-4;
        cfg.R *= 0.01;
        cfg.dt = 0.01;
        cfg.q0 = q0;
        return cfg;
    };

    // The truth is an EXACT fixed point. An instance seeded at the identity and
    // fed the gravity of the identity has an exactly zero innovation forever,
    // so it must not move at all -- bitwise, for the whole run. A correction
    // carrying any offset, or a propagation that drifts, fails this without a
    // threshold having to be chosen.
    {
        auto at_truth = build_mekf(configured(Eigen::Quaterniond::Identity()));
        const ctrlpp::Vector<double, 7> seeded = at_truth.state();
        for(int k = 0; k < 500; ++k)
        {
            at_truth.predict(ctrlpp::Vector<double, 3>::Zero());
            REQUIRE(at_truth.update(gravity_of_identity()).has_value());
            REQUIRE(at_truth.state() == seeded);
            REQUIRE(std::abs(at_truth.attitude().norm() - 1.0) <= unit_norm_budget);
        }
    }

    // The tilt keeps falling, across checkpoints an order of magnitude apart.
    // This is the statement the old fitted bound could not make: a filter that
    // converges part of the way and then stalls at a floor passes any single
    // threshold below its own starting error, and fails this ladder.
    auto filter = build_mekf(configured(
        Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitX()))));

    double at_100 = 0.0;
    double at_500 = 0.0;
    double at_1000 = 0.0;
    double at_2000 = 0.0;

    for(int k = 1; k <= 2000; ++k)
    {
        filter.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_of_identity()).has_value());
        REQUIRE(std::abs(filter.attitude().norm() - 1.0) <= unit_norm_budget);

        const double tilt = rotation_angle(filter.attitude());
        if(k == 100)
            at_100 = tilt;
        if(k == 500)
            at_500 = tilt;
        if(k == 1000)
            at_1000 = tilt;
        if(k == 2000)
            at_2000 = tilt;
    }

    CHECK(at_100 < initial_tilt);
    CHECK(at_500 < at_100);
    CHECK(at_1000 < at_500);
    CHECK(at_2000 < at_1000);

    // A gravity measurement is invariant under rotation about gravity, so an
    // initial tilt about x and one about y are the same problem in different
    // coordinates and the filter must produce the same tilt trajectory. It does
    // so BITWISE, which is a far sharper statement than any tolerance: an
    // axis-asymmetric measurement Jacobian, gain or reset Jacobian breaks it
    // immediately.
    auto about_y = build_mekf(configured(
        Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitY()))));
    auto about_x = build_mekf(configured(
        Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitX()))));

    for(int k = 0; k < 500; ++k)
    {
        about_x.predict(ctrlpp::Vector<double, 3>::Zero());
        about_y.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(about_x.update(gravity_of_identity()).has_value());
        REQUIRE(about_y.update(gravity_of_identity()).has_value());
        REQUIRE(rotation_angle(about_x.attitude()) == rotation_angle(about_y.attitude()));
    }
}

TEST_CASE("MEKF near-half-turn yaw error is unobservable from gravity",
          "[mekf][hardening][robustness]")
{
    // The fixture starts within a hundredth of a radian of a half turn ABOUT
    // THE GRAVITY AXIS. A gravity measurement is invariant under exactly that
    // rotation -- h(q0) is bitwise h(identity) -- so this is not a case about
    // the attitude branch cut at all: it is a case about an error the sensor
    // cannot see. Nothing here may assert convergence; the filter is correct
    // precisely because it does not move.
    const Eigen::Quaterniond q_init(Eigen::AngleAxisd(
        std::numbers::pi - 0.01, Eigen::Vector3d::UnitZ()));

    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    cfg.dt = 0.01;
    cfg.q0 = q_init;

    auto filter = build_mekf(cfg);

    // The premise the whole case rests on, asserted rather than assumed.
    CHECK(gravity_measurement{}(q_init.normalized(), ctrlpp::Vector<double, 3>::Zero())
          == gravity_of_identity());

    const ctrlpp::Vector<double, 7> seeded = filter.state();
    const double initial_yaw_variance = filter.covariance()(2, 2);
    double previous_yaw_variance = initial_yaw_variance;

    for(int k = 0; k < 200; ++k)
    {
        filter.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_of_identity()).has_value());

        // Exactly zero, not small: the predicted measurement is the same vector
        // as the supplied one, so their difference is a difference of identical
        // doubles. The correction is therefore exactly zero and the estimate is
        // BITWISE what it was seeded with, 200 steps later.
        REQUIRE(filter.innovation() == ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(filter.state() == seeded);
        REQUIRE(std::abs(filter.attitude().norm() - 1.0) <= unit_norm_budget);

        // And the filter says so: the variance of the direction it cannot
        // observe grows at every single step, because the propagation injects
        // process noise there and no update ever removes any.
        REQUIRE(filter.covariance()(2, 2) > previous_yaw_variance);
        previous_yaw_variance = filter.covariance()(2, 2);
    }

    // The two directions gravity DOES observe went the other way, ending far
    // below where they started. Stating both halves is what makes this a
    // statement about observability rather than about arithmetic surviving.
    CHECK(filter.covariance()(0, 0) < 1.0);
    CHECK(filter.covariance()(1, 1) < 1.0);
    CHECK(filter.covariance()(2, 2) > initial_yaw_variance);
    // The fixture is symmetric about the gravity axis, so the two observed
    // variances are equal bitwise.
    CHECK(filter.covariance()(0, 0) == filter.covariance()(1, 1));
}
