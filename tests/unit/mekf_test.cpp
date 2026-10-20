#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/lie/so3.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
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
        REQUIRE(filter.update(gravity_world).has_value());
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
        REQUIRE(filter.update(z).has_value());
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
    REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());

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
        REQUIRE(filter.update(gravity_world).has_value());
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
        REQUIRE(filter.update(gravity_world).has_value());
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

TEST_CASE("mekf create refuses a quaternion that is not unit within rounding", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{2.0, 0.0, 0.0, 0.0};

    auto filter = Mekf3::create(gravity_measurement{}, cfg);
    REQUIRE_FALSE(filter.has_value());
    REQUIRE(filter.error() == filter_error::non_unit_quaternion);
}

TEST_CASE("mekf create normalizes a quaternion whose norm drifts within rounding", "[mekf]")
{
    // Each component 0.5 + eps squares to 0.25 + eps once rounded, so the norm
    // is exactly 1 + 2 eps and the normalized quaternion exactly 0.5 each.
    const double component = 0.5 + std::numeric_limits<double>::epsilon();
    mekf_config<double, 3, 3> drifted;
    drifted.q0 = Eigen::Quaternion<double>{component, component, component, component};
    REQUIRE(drifted.q0.norm() == 1.0 + 2.0 * std::numeric_limits<double>::epsilon());
    mekf_config<double, 3, 3> unit = drifted;
    unit.q0 = drifted.q0.normalized();

    auto created = Mekf3::create(gravity_measurement{}, drifted);
    auto twin_created = Mekf3::create(gravity_measurement{}, unit);
    REQUIRE(created.has_value());
    REQUIRE(twin_created.has_value());
    auto& filter = *created;
    auto& twin = *twin_created;
    REQUIRE((filter.attitude().coeffs().array() == unit.q0.coeffs().array()).all());

    for(int k = 0; k < 10; ++k)
    {
        const Vector<double, 3> omega{0.2, -0.1, 0.3};
        const Vector<double, 3> z = Vector<double, 3>{0.05, 0.01 * k, 1.0}.normalized();
        filter.predict(omega);
        twin.predict(omega);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(twin.update(z).has_value());
    }
    REQUIRE((filter.attitude().coeffs().array() == twin.attitude().coeffs().array()).all());
    REQUIRE((filter.bias().array() == twin.bias().array()).all());
    REQUIRE((filter.covariance().array() == twin.covariance().array()).all());
}

TEST_CASE("mekf create refuses a noise, bias or initial covariance it cannot use", "[mekf]")
{
    const Matrix<double, 6, 6> indefinite = (Vector<double, 6>() << 1.0, 1.0, 1.0, 1.0, 1.0, -1.0).finished().asDiagonal();
    mekf_config<double, 3, 3> cfg;
    filter_error expected{};

    SECTION("a degenerate quaternion is reported ahead of the matrices")
    {
        cfg.q0 = Eigen::Quaternion<double>{0.0, 0.0, 0.0, 0.0};
        cfg.Q = indefinite;
        expected = filter_error::degenerate_quaternion;
    }

    SECTION("non-finite process noise")
    {
        cfg.Q(4, 4) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_process_noise;
    }

    SECTION("indefinite process noise")
    {
        cfg.Q = indefinite;
        expected = filter_error::indefinite_process_noise;
    }

    SECTION("non-finite measurement noise")
    {
        cfg.R(0, 1) = std::numeric_limits<double>::infinity();
        expected = filter_error::non_finite_measurement_noise;
    }

    SECTION("indefinite measurement noise")
    {
        cfg.R(2, 2) = -1.0;
        expected = filter_error::indefinite_measurement_noise;
    }

    SECTION("a non-unit quaternion is reported ahead of the matrices")
    {
        cfg.q0 = Eigen::Quaternion<double>{2.0, 0.0, 0.0, 0.0};
        cfg.Q = indefinite;
        expected = filter_error::non_unit_quaternion;
    }

    SECTION("a non-finite process noise is reported ahead of a non-finite initial bias")
    {
        cfg.Q(0, 0) = std::numeric_limits<double>::quiet_NaN();
        cfg.b0(1) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_process_noise;
    }

    SECTION("a non-finite measurement noise is reported ahead of a non-finite initial bias")
    {
        cfg.R(1, 1) = std::numeric_limits<double>::quiet_NaN();
        cfg.b0(1) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_measurement_noise;
    }

    SECTION("a non-finite initial bias is reported ahead of an indefinite initial covariance")
    {
        cfg.b0(1) = std::numeric_limits<double>::quiet_NaN();
        cfg.P0 = indefinite;
        expected = filter_error::non_finite_initial_bias;
    }

    SECTION("a non-finite initial bias is reported ahead of a non-finite initial covariance")
    {
        cfg.b0(2) = std::numeric_limits<double>::quiet_NaN();
        cfg.P0(0, 0) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_initial_bias;
    }

    SECTION("non-finite initial covariance")
    {
        cfg.P0(0, 0) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_initial_covariance;
    }

    SECTION("indefinite initial covariance")
    {
        cfg.P0 = indefinite;
        expected = filter_error::indefinite_initial_covariance;
    }

    SECTION("a degenerate quaternion is reported ahead of an asymmetric matrix")
    {
        cfg.q0 = Eigen::Quaternion<double>{0.0, 0.0, 0.0, 0.0};
        cfg.R(0, 2) = 0.5;
        expected = filter_error::degenerate_quaternion;
    }

    SECTION("structurally asymmetric process noise")
    {
        cfg.Q(5, 1) = 0.5;
        expected = filter_error::asymmetric_process_noise;
    }

    SECTION("structurally asymmetric measurement noise")
    {
        cfg.R(0, 2) = 0.5;
        expected = filter_error::asymmetric_measurement_noise;
    }

    SECTION("an asymmetric measurement noise is reported ahead of a non-finite initial bias")
    {
        cfg.R(0, 2) = 0.5;
        cfg.b0(0) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::asymmetric_measurement_noise;
    }

    SECTION("structurally asymmetric initial covariance")
    {
        cfg.P0(3, 4) = 0.5;
        expected = filter_error::asymmetric_initial_covariance;
    }

    const auto rejected = Mekf3::create(gravity_measurement{}, cfg);
    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == expected);
}

TEST_CASE("mekf create stores a drift-asymmetric matrix as its midpoint symmetrization", "[mekf]")
{
    mekf_config<double, 3, 3> drifted;
    drifted.Q = Matrix<double, 6, 6>::Identity() * 1e-4;
    drifted.Q(1, 4) = 1e-6;
    drifted.Q(4, 1) = std::nextafter(1e-6, 1.0);
    drifted.R = Matrix<double, 3, 3>::Identity() * 0.05;
    drifted.P0 = Matrix<double, 6, 6>::Identity() * 0.2;
    drifted.P0(0, 3) = 0.01;
    drifted.P0(3, 0) = std::nextafter(0.01, 0.0);

    mekf_config<double, 3, 3> symmetric = drifted;
    symmetric.Q(1, 4) = symmetric.Q(4, 1) = drifted.Q(1, 4) + (drifted.Q(4, 1) - drifted.Q(1, 4)) / 2.0;
    symmetric.P0(0, 3) = symmetric.P0(3, 0) = drifted.P0(0, 3) + (drifted.P0(3, 0) - drifted.P0(0, 3)) / 2.0;

    auto created = Mekf3::create(gravity_measurement{}, drifted);
    auto twin_created = Mekf3::create(gravity_measurement{}, symmetric);
    REQUIRE(created.has_value());
    REQUIRE(twin_created.has_value());
    auto& filter = *created;
    auto& twin = *twin_created;
    REQUIRE((filter.covariance().array() == symmetric.P0.array()).all());

    for(int k = 0; k < 10; ++k)
    {
        const Vector<double, 3> omega{0.2, -0.1, 0.3};
        const Vector<double, 3> z = Vector<double, 3>{0.05, 0.01 * k, 1.0}.normalized();
        filter.predict(omega);
        twin.predict(omega);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(twin.update(z).has_value());
    }
    REQUIRE((filter.attitude().coeffs().array() == twin.attitude().coeffs().array()).all());
    REQUIRE((filter.bias().array() == twin.bias().array()).all());
    REQUIRE((filter.covariance().array() == twin.covariance().array()).all());
}

TEST_CASE("mekf create accepts a singular positive semidefinite process noise", "[mekf]")
{
    mekf_config<double, 3, 3> cfg;
    cfg.Q.setZero();

    const auto accepted = Mekf3::create(gravity_measurement{}, cfg);
    REQUIRE(accepted.has_value());
    REQUIRE(accepted->covariance() == cfg.P0);
}
