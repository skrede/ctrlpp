#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/lie/so3.h"

#include "ctrlpp/estimation/sigma_points/so3_sigma_points.h"

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

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. The rejection
// cases below do not use this helper: they assert the specific enumerator.
auto make_filter(simple_rotation_dynamics dyn, const manifold_ukf_config<double, 3>& cfg) -> MukfType
{
    auto created = MukfType::create(std::move(dyn), gravity_meas{}, cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

constexpr double mukf_eps = std::numeric_limits<double>::epsilon();

struct tangent_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> Vector<double, 3>
    {
        return so3::log(q);
    }
};

// q <- exp((3b - b |log q|^2) e_z). With Merwe alpha = 1, beta = 0, kappa = -1
// at n = 3 the center weight is -1/2, the six outer weights are 1/4 and the
// outer points of P = I sit at |log q|^2 = 2, so the center maps to the angle 3b
// and the outer points to b: the weighted mean angle is -3b/2 + 6b/4 = 0 and the
// weighted variance is -9b^2/2 + 6b^2/4 = -3b^2, from a positive definite prior.
struct z_spread_dynamics
{
    double b;

    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& /*omega*/) const -> Eigen::Quaternion<double>
    {
        const double angle = 3.0 * b - b * so3::log(q).squaredNorm();
        return so3::exp(Vector<double, 3>{0.0, 0.0, angle});
    }
};

using tangent_filter_t = manifold_ukf<double, 3, z_spread_dynamics, tangent_measurement>;

constexpr double z_spread_b = 0.1;

auto z_spread_filter(manifold_ukf_config<double, 3> cfg) -> tangent_filter_t
{
    auto strategy = so3_merwe_sigma_points<double>::try_create({.alpha = 1.0, .beta = 0.0, .kappa = -1.0});
    REQUIRE(strategy.has_value());
    cfg.Q.setZero();
    auto created = tangent_filter_t::create(z_spread_dynamics{z_spread_b}, tangent_measurement{}, cfg, *strategy);
    REQUIRE(created.has_value());
    return *std::move(created);
}

/// @brief The floor a repaired covariance is held to: the symmetric
/// eigensolver's backward error 3 eps max|P_ij|, the same scale the repair
/// clamps its spectrum to.
auto repair_floor(const Matrix<double, 3, 3>& P) -> double
{
    return 3.0 * mukf_eps * P.cwiseAbs().maxCoeff();
}

auto min_eigenvalue(const Matrix<double, 3, 3>& P) -> double
{
    return Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3>>(P).eigenvalues().minCoeff();
}

} // namespace

TEST_CASE("manifold_ukf predict with zero angular velocity preserves attitude", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

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
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_world).has_value());
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
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

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

    Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("manifold_ukf geodesic mean converges", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 30; ++i)
    {
        filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
        REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());
    }

    auto q = filter.attitude();
    double qnorm = std::sqrt(q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
    REQUIRE_THAT(qnorm, WithinAbs(1.0, 1e-10));
}

TEST_CASE("manifold_ukf innovation is finite", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
    REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());

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
    auto filter = make_filter(dyn, cfg);

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
        REQUIRE(filter.update(z).has_value());
    }

    auto q_est = filter.attitude();
    Vector<double, 3> phi = so3::log(q_est);
    REQUIRE(std::abs(phi(2)) > 0.01); // Some rotation accumulated
}

// ---------------------------------------------------------------------------
// Edge / failure-path tests
// ---------------------------------------------------------------------------

TEST_CASE("manifold_ukf predict-only grows covariance", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q = Matrix<double, 3, 3>::Identity() * 0.01;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    auto P_before = filter.covariance();
    double trace_before = P_before.trace();

    for(int i = 0; i < 50; ++i)
        filter.predict(Vector<double, 3>::Zero());

    double trace_after = filter.covariance().trace();
    // Covariance should grow when only predicting (no measurement updates)
    REQUIRE(trace_after > trace_before);
}

TEST_CASE("manifold_ukf predict-only covariance stays symmetric PSD", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q = Matrix<double, 3, 3>::Identity() * 0.001;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 100; ++i)
        filter.predict(Vector<double, 3>{0.1, -0.05, 0.02});

    auto P = filter.covariance();
    double asym = (P - P.transpose()).norm();
    REQUIRE(asym < 1e-10);

    Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("manifold_ukf large initial error converges with updates", "[manifold_ukf]")
{
    // Start with ~170 degree error (near antipodal)
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = so3::exp(Vector<double, 3>{2.9, 0.0, 0.0}); // ~166 degrees
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.P0 *= 100.0; // Very uncertain
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 200; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_world).has_value());
    }

    auto q = filter.attitude();
    // Should have made progress toward identity (gravity-aligned)
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    // With large initial error the filter should at least reduce the error significantly
    REQUIRE(angle < 2.0);
}

TEST_CASE("manifold_ukf with very small process noise", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q = Matrix<double, 3, 3>::Identity() * 1e-12;
    cfg.R *= 0.1;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    // The filter should still work without numerical issues
    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());
    }

    auto q = filter.attitude();
    REQUIRE(std::isfinite(q.w()));
    REQUIRE(std::isfinite(q.x()));
    REQUIRE(std::isfinite(q.y()));
    REQUIRE(std::isfinite(q.z()));

    auto P = filter.covariance();
    REQUIRE(std::isfinite(P.norm()));
}

TEST_CASE("manifold_ukf with very large measurement noise trusts prediction", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R = Matrix<double, 3, 3>::Identity() * 1e6; // Huge measurement noise
    cfg.P0 *= 0.001; // Very confident initial state
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    // Give a measurement that disagrees with identity orientation
    for(int i = 0; i < 20; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        // Wrong direction, but a perfectly finite measurement, so the step
        // succeeds and the filter is expected to disagree with it.
        REQUIRE(filter.update(Vector<double, 3>{0.5, 0.5, 0.5}).has_value());
    }

    // Filter should mostly ignore measurement due to high R
    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.5); // Still close to identity
}

TEST_CASE("manifold_ukf geodesic mean with max_iter=1", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.geodesic_mean_max_iter = 1;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    // Should still produce finite results even with single geodesic iteration
    for(int i = 0; i < 20; ++i)
    {
        filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
        REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());
    }

    auto q = filter.attitude();
    REQUIRE(std::isfinite(q.w()));
    double qnorm = std::sqrt(q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
    REQUIRE_THAT(qnorm, WithinAbs(1.0, 1e-10));
}

TEST_CASE("manifold_ukf with loose geodesic tolerance", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.geodesic_mean_tol = 1.0; // Very loose -- should converge immediately
    cfg.Q *= 1e-4;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 10; ++i)
        filter.predict(Vector<double, 3>{0.01, -0.01, 0.02});

    auto q = filter.attitude();
    REQUIRE(std::isfinite(q.w()));
}

TEST_CASE("manifold_ukf state cache matches attitude quaternion", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    filter.predict(Vector<double, 3>{0.05, -0.03, 0.01});
    REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());

    auto q = filter.attitude();
    auto s = filter.state();

    // state_cache should be so3::to_vec(q) = [w, x, y, z]
    REQUIRE_THAT(s(0), WithinAbs(q.w(), 1e-15));
    REQUIRE_THAT(s(1), WithinAbs(q.x(), 1e-15));
    REQUIRE_THAT(s(2), WithinAbs(q.y(), 1e-15));
    REQUIRE_THAT(s(3), WithinAbs(q.z(), 1e-15));
}

TEST_CASE("manifold_ukf innovation decreases as filter converges", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-8;
    cfg.R *= 0.01;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> gravity{0.0, 0.0, 1.0};

    // Run several predict/update cycles to let the filter converge
    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity).has_value());
    }

    // After convergence, innovation should be small
    double innov_norm = filter.innovation().norm();
    REQUIRE(innov_norm < 0.1);
    REQUIRE(std::isfinite(innov_norm));
}

TEST_CASE("manifold_ukf create rejects a zero initial quaternion", "[manifold_ukf]")
{
    // Before the fallible factory existed, the constructor normalized the zero
    // quaternion directly and silently produced an all-NaN filter state.
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{0.0, 0.0, 0.0, 0.0};

    auto filter = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE_FALSE(filter.has_value());
    REQUIRE(filter.error() == filter_error::degenerate_quaternion);
}

TEST_CASE("manifold_ukf create refuses a quaternion that is not unit within rounding", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{2.0, 0.0, 0.0, 0.0};

    auto filter = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE_FALSE(filter.has_value());
    REQUIRE(filter.error() == filter_error::non_unit_quaternion);
}

TEST_CASE("manifold_ukf create normalizes a quaternion whose norm drifts within rounding", "[manifold_ukf]")
{
    // Each component 0.5 + eps squares to 0.25 + eps once rounded, so the norm
    // is exactly 1 + 2 eps and the normalized quaternion exactly 0.5 each.
    const double component = 0.5 + std::numeric_limits<double>::epsilon();
    manifold_ukf_config<double, 3> drifted;
    drifted.q0 = Eigen::Quaternion<double>{component, component, component, component};
    REQUIRE(drifted.q0.norm() == 1.0 + 2.0 * std::numeric_limits<double>::epsilon());
    manifold_ukf_config<double, 3> unit = drifted;
    unit.q0 = drifted.q0.normalized();

    auto created = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, drifted);
    auto twin_created = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, unit);
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
    REQUIRE((filter.covariance().array() == twin.covariance().array()).all());
}

TEST_CASE("manifold_ukf create refuses a noise or initial covariance that is not a covariance", "[manifold_ukf]")
{
    const Matrix<double, 3, 3> indefinite = Vector<double, 3>{1.0, 1.0, -1.0}.asDiagonal();
    manifold_ukf_config<double, 3> cfg;
    filter_error expected{};

    SECTION("a degenerate quaternion is reported ahead of the matrices")
    {
        cfg.q0 = Eigen::Quaternion<double>{0.0, 0.0, 0.0, 0.0};
        cfg.Q = indefinite;
        expected = filter_error::degenerate_quaternion;
    }

    SECTION("a non-unit quaternion is reported ahead of the matrices")
    {
        cfg.q0 = Eigen::Quaternion<double>{2.0, 0.0, 0.0, 0.0};
        cfg.Q = indefinite;
        expected = filter_error::non_unit_quaternion;
    }

    SECTION("non-finite process noise")
    {
        cfg.Q(0, 0) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_process_noise;
    }

    SECTION("indefinite process noise")
    {
        cfg.Q = indefinite;
        expected = filter_error::indefinite_process_noise;
    }

    SECTION("non-finite measurement noise")
    {
        cfg.R(1, 1) = std::numeric_limits<double>::quiet_NaN();
        expected = filter_error::non_finite_measurement_noise;
    }

    SECTION("indefinite measurement noise")
    {
        cfg.R = indefinite;
        expected = filter_error::indefinite_measurement_noise;
    }

    SECTION("non-finite initial covariance")
    {
        cfg.P0(2, 2) = std::numeric_limits<double>::infinity();
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
        cfg.Q(0, 1) = 0.5;
        expected = filter_error::degenerate_quaternion;
    }

    SECTION("structurally asymmetric process noise")
    {
        cfg.Q(0, 1) = 0.5;
        expected = filter_error::asymmetric_process_noise;
    }

    SECTION("structurally asymmetric measurement noise")
    {
        cfg.R(2, 0) = 0.5;
        expected = filter_error::asymmetric_measurement_noise;
    }

    SECTION("structurally asymmetric initial covariance")
    {
        cfg.P0(1, 2) = 0.5;
        expected = filter_error::asymmetric_initial_covariance;
    }

    const auto rejected = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == expected);
}

TEST_CASE("manifold_ukf create stores a drift-asymmetric matrix as its midpoint symmetrization", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> drifted;
    drifted.Q = Vector<double, 3>{1e-4, 2e-4, 1e-4}.asDiagonal();
    drifted.Q(0, 2) = 2e-5;
    drifted.Q(2, 0) = std::nextafter(2e-5, 1.0);
    drifted.R = Matrix<double, 3, 3>::Identity() * 0.05;
    drifted.P0 = Matrix<double, 3, 3>::Identity() * 0.3;
    drifted.P0(0, 1) = 0.05;
    drifted.P0(1, 0) = std::nextafter(0.05, 0.0);

    manifold_ukf_config<double, 3> symmetric = drifted;
    symmetric.Q(0, 2) = symmetric.Q(2, 0) = drifted.Q(0, 2) + (drifted.Q(2, 0) - drifted.Q(0, 2)) / 2.0;
    symmetric.P0(0, 1) = symmetric.P0(1, 0) = drifted.P0(0, 1) + (drifted.P0(1, 0) - drifted.P0(0, 1)) / 2.0;

    auto filter = make_filter(simple_rotation_dynamics{}, drifted);
    auto twin = make_filter(simple_rotation_dynamics{}, symmetric);
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
    REQUIRE((filter.covariance().array() == twin.covariance().array()).all());
}

TEST_CASE("manifold_ukf create accepts a singular positive semidefinite configuration", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q.setZero();
    cfg.P0 = Vector<double, 3>{1.0, 1.0, 0.0}.asDiagonal();

    const auto accepted = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE(accepted.has_value());
    REQUIRE(accepted->covariance() == cfg.P0);
}

TEST_CASE("manifold_ukf refuses an indefinite initial covariance at construction", "[manifold_ukf][repair]")
{
    // The repaired-update property is covered by the prediction-route case below.
    manifold_ukf_config<double, 3> cfg;
    cfg.Q.setZero();
    cfg.P0 = Vector<double, 3>{1.0, 1.0, -0.5}.asDiagonal();
    const auto rejected = manifold_ukf<double, 3, simple_rotation_dynamics, tangent_measurement>::create(simple_rotation_dynamics{}, tangent_measurement{}, cfg);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == filter_error::indefinite_initial_covariance);
}

TEST_CASE("manifold_ukf update after a prediction-driven repair carries the repaired covariance", "[manifold_ukf][repair]")
{
    auto filter = z_spread_filter(manifold_ukf_config<double, 3>{});

    filter.predict(Vector<double, 3>::Zero());

    // Seven weighted terms, each through one log, one product and one sum, at
    // most the center's angle 3b for the mean and its square for the variance.
    constexpr double sigma_ops = 3.0 * 7.0;
    const double angle_scale = 3.0 * z_spread_b;
    CAPTURE(filter.attitude().coeffs(), filter.covariance());
    REQUIRE(so3::log(filter.attitude()).norm() <= sigma_ops * mukf_eps * angle_scale);
    REQUIRE_THAT(filter.covariance()(2, 2), WithinAbs(-3.0 * z_spread_b * z_spread_b, sigma_ops * mukf_eps * angle_scale * angle_scale));
    REQUIRE(filter.health() == manifold_ukf_health::covariance_repaired);

    REQUIRE(filter.update(Vector<double, 3>::Zero()).has_value());

    CAPTURE(filter.covariance());
    REQUIRE(min_eigenvalue(filter.covariance()) >= -repair_floor(filter.covariance()));
    REQUIRE(filter.health() == manifold_ukf_health::covariance_repaired);
}

namespace
{

// Non-finite for every sigma point rotated off the z axis, so an update whose
// sigma set spreads about x or y produces a non-finite result.
struct off_axis_non_finite_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> Vector<double, 3>
    {
        const Vector<double, 3> phi = so3::log(q);
        if(phi.x() != 0.0 || phi.y() != 0.0)
            return Vector<double, 3>::Constant(std::numeric_limits<double>::quiet_NaN());
        return phi;
    }
};

struct non_finite_rotation_dynamics
{
    auto operator()(const Eigen::Quaternion<double>& /*q*/, const Vector<double, 3>& /*omega*/) const -> Eigen::Quaternion<double>
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return Eigen::Quaternion<double>{nan, nan, nan, nan};
    }
};

}

TEST_CASE("manifold_ukf predict reports a non-finite prediction before the next update", "[manifold_ukf][repair]")
{
    using filter_t = manifold_ukf<double, 3, non_finite_rotation_dynamics, gravity_meas>;
    auto created = filter_t::create(non_finite_rotation_dynamics{}, gravity_meas{}, manifold_ukf_config<double, 3>{});
    REQUIRE(created.has_value());
    auto filter = *std::move(created);

    filter.predict(Vector<double, 3>::Zero());
    REQUIRE(filter.health() == manifold_ukf_health::non_finite_estimate);

    const auto refused = filter.update(Vector<double, 3>::UnitZ());
    REQUIRE_FALSE(refused.has_value());
    REQUIRE(refused.error() == manifold_ukf_update_error::non_finite_state);
    REQUIRE(filter.health() == manifold_ukf_health::non_finite_estimate);
}

TEST_CASE("manifold_ukf update rolled back after an indefinite prediction keeps the rung predict latched", "[manifold_ukf][repair]")
{
    using filter_t = manifold_ukf<double, 3, z_spread_dynamics, off_axis_non_finite_measurement>;
    auto strategy = so3_merwe_sigma_points<double>::try_create({.alpha = 1.0, .beta = 0.0, .kappa = -1.0});
    REQUIRE(strategy.has_value());
    manifold_ukf_config<double, 3> cfg;
    cfg.Q.setZero();
    auto created = filter_t::create(z_spread_dynamics{z_spread_b}, off_axis_non_finite_measurement{}, cfg, *strategy);
    REQUIRE(created.has_value());
    auto filter = *std::move(created);

    filter.predict(Vector<double, 3>::Zero());
    REQUIRE(filter.health() == manifold_ukf_health::covariance_repaired);
    const auto attitude_before = filter.attitude().coeffs().eval();
    const auto covariance_before = filter.covariance();

    const auto refused = filter.update(Vector<double, 3>::Zero());
    REQUIRE_FALSE(refused.has_value());
    REQUIRE(refused.error() == manifold_ukf_update_error::non_finite_result);
    CHECK(filter.attitude().coeffs() == attitude_before);
    CHECK(filter.covariance() == covariance_before);
    REQUIRE(filter.health() == manifold_ukf_health::covariance_repaired);
}

TEST_CASE("manifold_ukf health latches a repair and a non-converged mean by severity", "[manifold_ukf][repair]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.geodesic_mean_max_iter = 1;

    SECTION("a later repair does not lower a non-converged mean")
    {
        auto filter = z_spread_filter(cfg);
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.health() == manifold_ukf_health::mean_not_converged);

        REQUIRE(filter.update(Vector<double, 3>::Zero()).has_value());
        REQUIRE(min_eigenvalue(filter.covariance()) >= -repair_floor(filter.covariance()));
        REQUIRE(filter.health() == manifold_ukf_health::mean_not_converged);
    }

    // The unpivoted Cholesky factorization fails on this singular prior, so the
    // first update repairs it.
    SECTION("a later non-converged mean raises a repair")
    {
        cfg.P0 = Vector<double, 3>{1.0, 1.0, 0.0}.asDiagonal();
        auto filter = z_spread_filter(cfg);
        REQUIRE(filter.update(Vector<double, 3>::Zero()).has_value());
        REQUIRE(filter.health() == manifold_ukf_health::covariance_repaired);

        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.health() == manifold_ukf_health::mean_not_converged);
    }
}
