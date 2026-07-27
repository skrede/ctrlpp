#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/lie/so3.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <utility>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

struct gravity_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& /*b*/) const -> Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

struct analytical_gravity_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& /*b*/) const -> Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }

    auto jacobian(const Eigen::Quaternion<double>& q, const Vector<double, 3>& /*b*/) const -> Matrix<double, 3, 6>
    {
        // Numerical approximation of the Jacobian via finite differences
        constexpr double eps = 1e-6;
        Matrix<double, 3, 6> H;
        Vector<double, 3> z0 = q.toRotationMatrix().transpose().col(2);

        for(int j = 0; j < 3; ++j)
        {
            Vector<double, 3> delta = Vector<double, 3>::Zero();
            delta(j) = eps;
            Eigen::Quaternion<double> qp = (q * so3::exp(delta)).normalized();
            Vector<double, 3> zp = qp.toRotationMatrix().transpose().col(2);
            H.col(j) = (zp - z0) / eps;
        }
        // Bias columns are zero (measurement does not depend on bias)
        H.rightCols<3>().setZero();
        return H;
    }
};

using Mekf3 = mekf<double, 3, 3, gravity_measurement>;
using MekfAnalytical = mekf<double, 3, 3, analytical_gravity_measurement>;

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. The rejection
// cases below do not use this helper: they assert the specific enumerator.
template <typename Filter, typename Measurement, typename Config>
auto make_filter(Measurement measurement, const Config& cfg) -> Filter
{
    auto created = Filter::create(std::move(measurement), cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

} // namespace

TEST_CASE("mekf predict with zero angular velocity preserves attitude", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-6;
    auto filter = make_filter<Mekf3>(gravity_measurement{}, cfg);

    for(int i = 0; i < 100; ++i)
        filter.predict(Vector<double, 3>::Zero());

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.01);
}

TEST_CASE("mekf predict integrates constant angular velocity", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.dt = 0.01;
    cfg.Q *= 1e-8;
    auto filter = make_filter<Mekf3>(gravity_measurement{}, cfg);

    Vector<double, 3> omega;
    omega << 0.0, 0.0, 0.1;

    for(int i = 0; i < 100; ++i)
        filter.predict(omega);

    // Expected total rotation: 100 * 0.01 * 0.1 = 0.1 rad about z
    auto q = filter.attitude();
    Vector<double, 3> phi = so3::log(q);
    REQUIRE_THAT(phi(2), WithinAbs(0.1, 0.05));
}

TEST_CASE("mekf update corrects attitude toward measurement", "[mekf]")
{
    // Start with 10-degree error about x-axis
    mekf_config<double, 3, 3> cfg;
    cfg.q0 = so3::exp(Vector<double, 3>{0.174, 0.0, 0.0});
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    auto filter = make_filter<Mekf3>(gravity_measurement{}, cfg);

    Vector<double, 3> gravity_world;
    gravity_world << 0.0, 0.0, 1.0;

    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        filter.update(gravity_world);
    }

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.05);
}

TEST_CASE("mekf covariance stays symmetric and PSD", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    auto filter = make_filter<Mekf3>(gravity_measurement{}, cfg);

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

    Eigen::SelfAdjointEigenSolver<Matrix<double, 6, 6>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("mekf innovation is finite after update", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    auto filter = make_filter<Mekf3>(gravity_measurement{}, cfg);

    filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
    filter.update(Vector<double, 3>{0.0, 0.0, 1.0});

    REQUIRE(std::isfinite(filter.innovation().norm()));
}

TEST_CASE("mekf bias estimation converges", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.dt = 0.01;
    cfg.Q = Matrix<double, 6, 6>::Identity() * 1e-4;
    cfg.Q.bottomRightCorner<3, 3>() *= 10.0; // More process noise on bias
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    auto filter = make_filter<Mekf3>(gravity_measurement{}, cfg);

    Vector<double, 3> true_bias{0.01, -0.02, 0.005};
    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 500; ++i)
    {
        // Gyro reads zero angular velocity + bias
        filter.predict(true_bias);
        filter.update(gravity_world);
    }

    auto est_bias = filter.bias();
    REQUIRE((est_bias - true_bias).norm() < 0.05);
}

TEST_CASE("mekf satisfies ObserverPolicy concept", "[mekf]")
{
    static_assert(ObserverPolicy<Mekf3>);
    static_assert(CovarianceObserver<Mekf3>);
}

TEST_CASE("mekf with analytical Jacobian measurement model", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.q0 = so3::exp(Vector<double, 3>{0.174, 0.0, 0.0});
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    auto filter = make_filter<MekfAnalytical>(analytical_gravity_measurement{}, cfg);

    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        filter.update(gravity_world);
    }

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.05);
}

TEST_CASE("mekf create rejects a zero initial quaternion", "[mekf]")
{
    // Before the fallible factory existed, the constructor normalized the zero
    // quaternion directly and silently produced an all-NaN filter state.
    mekf_config<double, 3, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{0.0, 0.0, 0.0, 0.0};

    auto filter = Mekf3::create(gravity_measurement{}, cfg);
    REQUIRE_FALSE(filter.has_value());
    REQUIRE(filter.error() == filter_error::degenerate_quaternion);
}

TEST_CASE("mekf create accepts a non-unit nonzero quaternion", "[mekf]")
{
    // The guard rejects only degeneracy: any finite nonzero quaternion is
    // normalized onto the unit sphere at construction.
    mekf_config<double, 3, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{2.0, 0.0, 0.0, 0.0};

    auto filter = Mekf3::create(gravity_measurement{}, cfg);
    REQUIRE(filter.has_value());
    REQUIRE((filter->attitude().coeffs().array() == cfg.q0.normalized().coeffs().array()).all());
}
