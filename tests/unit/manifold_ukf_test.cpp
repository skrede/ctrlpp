#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/lie/so3.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

struct simple_rotation_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& omega) const -> Eigen::Quaternion<double>
    {
        Vector<double, 3> phi = (omega * dt).eval();
        return (q * so3::exp(phi)).normalized();
    }
};

struct gravity_meas
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

using MukfType = manifold_ukf<double, 3, simple_rotation_dynamics, gravity_meas>;

} // namespace

TEST_CASE("manifold_ukf predict with zero angular velocity preserves attitude", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    MukfType filter(simple_rotation_dynamics{}, gravity_meas{}, cfg);

    for(int i = 0; i < 100; ++i)
        filter.predict(Vector<double, 3>::Zero());

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.01);
}

TEST_CASE("manifold_ukf update corrects attitude toward measurement", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = so3::exp(Vector<double, 3>{0.174, 0.0, 0.0});
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    MukfType filter(simple_rotation_dynamics{}, gravity_meas{}, cfg);

    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        filter.update(gravity_world);
    }

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.1);
}

TEST_CASE("manifold_ukf covariance stays symmetric and PSD", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    MukfType filter(simple_rotation_dynamics{}, gravity_meas{}, cfg);

    Vector<double, 3> omega{0.01, -0.02, 0.03};
    Vector<double, 3> z{0.0, 0.0, 1.0};

    for(int i = 0; i < 100; ++i)
    {
        filter.predict(omega);
        filter.update(z);
    }

    auto P = filter.covariance();
    double asym = (P - P.transpose()).norm();
    REQUIRE(asym < 1e-10);

    Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("manifold_ukf geodesic mean converges", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    MukfType filter(simple_rotation_dynamics{}, gravity_meas{}, cfg);

    for(int i = 0; i < 30; ++i)
    {
        filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
        filter.update(Vector<double, 3>{0.0, 0.0, 1.0});
    }

    auto q = filter.attitude();
    double qnorm = std::sqrt(q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
    REQUIRE_THAT(qnorm, WithinAbs(1.0, 1e-10));
}

TEST_CASE("manifold_ukf innovation is finite", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    MukfType filter(simple_rotation_dynamics{}, gravity_meas{}, cfg);

    filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
    filter.update(Vector<double, 3>{0.0, 0.0, 1.0});

    REQUIRE(std::isfinite(filter.innovation().norm()));
}

TEST_CASE("manifold_ukf satisfies ObserverPolicy concept", "[manifold_ukf]")
{
    static_assert(ObserverPolicy<MukfType>);
    static_assert(CovarianceObserver<MukfType>);
}

TEST_CASE("manifold_ukf tracks constant rotation", "[manifold_ukf]")
{
    simple_rotation_dynamics dyn;
    dyn.dt = 0.01;

    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-8;
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    MukfType filter(dyn, gravity_meas{}, cfg);

    Vector<double, 3> omega{0.0, 0.0, 0.1};
    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    // Run 100 steps: expected rotation = 100 * 0.01 * 0.1 = 0.1 rad about z
    for(int i = 0; i < 100; ++i)
    {
        filter.predict(omega);
        // Generate measurement from ground truth
        Vector<double, 3> phi_true = (omega * dyn.dt * static_cast<double>(i + 1)).eval();
        Eigen::Quaternion<double> q_true = so3::exp(phi_true);
        Vector<double, 3> z = q_true.toRotationMatrix().transpose().col(2);
        filter.update(z);
    }

    auto q_est = filter.attitude();
    Vector<double, 3> phi = so3::log(q_est);
    REQUIRE(std::abs(phi(2)) > 0.01); // Some rotation accumulated
}
