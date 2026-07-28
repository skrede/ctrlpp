#include "hardening_helpers.h"
#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace {

struct simple_rotation_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double>& q,
                    const ctrlpp::Vector<double, 3>& omega) const -> Eigen::Quaternion<double>
    {
        ctrlpp::Vector<double, 3> phi = (omega * dt).eval();
        return (q * ctrlpp::so3::exp(phi)).normalized();
    }
};

struct gravity_meas
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> ctrlpp::Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

using MukfType = ctrlpp::manifold_ukf<double, 3, simple_rotation_dynamics, gravity_meas>;

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here.
auto make_filter(const ctrlpp::manifold_ukf_config<double, 3>& cfg) -> MukfType
{
    auto created = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

auto quat_angle(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) -> double
{
    return 2.0 * std::acos(std::min(1.0, std::abs(q1.dot(q2))));
}

}

TEST_CASE("Manifold UKF non-unit quaternion input normalizes",
          "[manifold_ukf][hardening][negative]")
{
    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;

    // Non-unit quaternion as initial condition
    Eigen::Quaterniond q_init(2.0, 0.0, 0.0, 0.0);
    cfg.q0 = q_init;

    auto filter = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(omega);

    ctrlpp::Vector<double, 3> z;
    z << 0.0, 0.0, 1.0;
    REQUIRE(filter.update(z).has_value());

    CHECK(std::isfinite(filter.state()[0]));
}

TEST_CASE("Manifold UKF NaN rotation measurement is rejected without touching the estimate",
          "[manifold_ukf][hardening][negative]")
{
    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;

    auto filter = make_filter(cfg);
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(omega);
    reference.predict(omega);

    // Snapshot immediately before the poisoned step.
    const ctrlpp::Vector<double, 4> state_before = filter.state();
    const ctrlpp::Matrix<double, 3, 3> P_before = filter.covariance();

    ctrlpp::Vector<double, 3> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0;

    const auto rejected = filter.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::manifold_ukf_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(filter.state() == state_before);
    // The covariance was already measurement-independent before the guard
    // existed -- the tangent reduction P - K*S*K^T and its reset Jacobian are
    // built from the sigma points, the gain and the correction, never from z --
    // so this half of the invariant is structural. The state half is what the
    // guard adds.
    CHECK(filter.covariance() == P_before);
    // On a manifold the attitude quaternion's unit norm is an extra invariant,
    // and the rejection preserves it exactly rather than approximately.
    CHECK(filter.attitude().norm() == 1.0);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(filter.health() == ctrlpp::manifold_ukf_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    ctrlpp::Vector<double, 3> z_good;
    z_good << 0.0, 0.0, 1.0;
    REQUIRE(filter.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(filter.state() == reference.state());
    CHECK(filter.covariance() == reference.covariance());
}

TEST_CASE("Manifold UKF covariance stays PD over 1000 steps",
          "[manifold_ukf][hardening][stability]")
{
    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;

    auto filter = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();
    bool all_pd = true;

    for(int k = 0; k < 1000; ++k)
    {
        filter.predict(omega);
        ctrlpp::Vector<double, 3> z;
        z << 0.0, 0.0, 1.0;
        REQUIRE(filter.update(z).has_value());

        Eigen::SelfAdjointEigenSolver<ctrlpp::Matrix<double, 3, 3>> eigsolver(filter.covariance());
        for(int i = 0; i < 3; ++i)
        {
            if(eigsolver.eigenvalues()(i) < -1e-10)
            {
                all_pd = false;
                break;
            }
        }
        if(!all_pd)
            break;
    }

    REQUIRE(all_pd);
}

TEST_CASE("Manifold UKF attitude converges for slow rotation",
          "[manifold_ukf][hardening][convergence]")
{
    // Start tilted
    auto q_init = Eigen::Quaterniond(Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitY()));

    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.q0 = q_init;

    auto filter = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();

    for(int k = 0; k < 500; ++k)
    {
        filter.predict(omega);
        // Gravity measurement for identity attitude
        ctrlpp::Vector<double, 3> z;
        z << 0.0, 0.0, 1.0;
        REQUIRE(filter.update(z).has_value());
    }

    auto q_est = filter.attitude();
    double angle_error = quat_angle(q_est, Eigen::Quaterniond::Identity());
    REQUIRE(angle_error < 0.1);
}

TEST_CASE("Manifold UKF extreme rotation near gimbal lock",
          "[manifold_ukf][hardening][robustness]")
{
    // Start near 180 degrees
    auto q_init = Eigen::Quaterniond(Eigen::AngleAxisd(
        std::numbers::pi - 0.05, Eigen::Vector3d::UnitX()));

    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    cfg.q0 = q_init;

    auto filter = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();
    bool all_finite = true;

    for(int k = 0; k < 200; ++k)
    {
        filter.predict(omega);
        ctrlpp::Vector<double, 3> z;
        z << 0.0, 0.0, 1.0;
        REQUIRE(filter.update(z).has_value());

        if(!std::isfinite(filter.state()[0]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}

TEST_CASE("Manifold UKF near pi rotation triggers hemisphere flip",
          "[manifold_ukf][hardening][coverage]")
{
    // Initial quaternion near 180 degrees around z-axis -- sigma points will
    // span both hemispheres, exercising the qi.dot(q_mean) < 0 branch.
    Eigen::Quaterniond q_near_pi(Eigen::AngleAxisd(3.0, Eigen::Vector3d::UnitZ()));

    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.q0 = q_near_pi;
    cfg.P0 *= 0.5; // Large covariance -> wide sigma spread
    cfg.Q *= 0.01;
    cfg.R *= 0.1;

    auto filter = make_filter(cfg);

    Eigen::Vector3d omega(0.1, 0.0, 0.5);

    for(int k = 0; k < 50; ++k)
    {
        filter.predict(omega);
        Eigen::Vector3d z = q_near_pi.toRotationMatrix().transpose().col(2);
        REQUIRE(filter.update(z).has_value());
    }

    // Filter should remain finite through hemisphere-spanning sigma points
    REQUIRE(std::isfinite(filter.state()[0]));
}
